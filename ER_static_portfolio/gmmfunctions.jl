# ==========================================================================
# 4D MODEL: GMM/SMM estimation functions
# ==========================================================================

using DelimitedFiles
using Printf
using LinearAlgebra

# ==========================================================================
# Build full parameter vector from estimated subset
# ==========================================================================
function buildparam(p::Vector{Float64})
    # p contains the estimated parameters: [nu, F^d, kappa, chi, F^t]
    # Returns full 17-element parameter vector
    
    pea = zeros(sz.nop)
    
    # Fixed parameters (calibrated)
    pea[1] = 0.98       # beta (discount factor)
    pea[2] = 0.025      # delta (depreciation)
    pea[3] = 0.9        # rho_e (ER persistence)
    pea[4] = 0.15       # sigma_e (ER volatility)
    pea[6] = 2.0        # gamma (risk aversion)
    pea[8] = 1.0        # wage
    pea[9] = 0.02       # r_star (dollar rate)
    pea[10] = 1.0       # pd (durable price)
    pea[12] = 0.25      # tau (tax rate)
    pea[13] = 1.0       # h (labor supply)
    pea[14] = 0.9       # rho_y (income persistence)
    pea[15] = 0.1       # sigma_y (income volatility)
    
    # Estimated parameters
    pea[5] = p[1]       # nu (non-durable share)
    pea[7] = p[2]       # F^d (durable fixed cost)
    pea[11] = p[3]      # kappa (dollar transaction cost)
    pea[16] = p[4]      # chi (maintenance effectiveness)
    pea[17] = p[5]      # F^t (time cost)
    
    return pea::Vector{Float64}
end


# ==========================================================================
# Main fcn: returns moments given parameters
# ==========================================================================
function fcn(p::Vector{Float64})
    pea = buildparam(p)
    moms = momentgen(pea)
    moms = moms[sz.pick]
    return moms::Vector{Float64}
end


# ==========================================================================
# fcn with GMM objective value
# ==========================================================================
function fcn(p::Vector{Float64}, fopt::Float64)
    pea = buildparam(p)
    
    simmoms = collect(momentgen(pea))
    
    if simmoms[1] > -99.0  # Valid simulation
        datamoms = vec(collect(readdlm(kst.MOMS_FILE)))
        sw = collect(readdlm(kst.W_FILE))
        momname = vec(collect(readdlm(kst.MNAME_FILE)))
        
        ch = sz.pick
        sw = sw[ch, ch]
        
        if settings.complicated
            w = I(size(sw, 1))
        else
            w = inv(sw)
        end
        
        datamoms = datamoms[ch]
        simmoms = simmoms[ch]
        momdiff = datamoms .- simmoms
        momname = momname[ch]
        
        bigQ = (momdiff' * w * momdiff)[1]
        
        # Save progress if improved
        if bigQ < fopt
            io = open(kst.PROGRESS_FILE, "w")
            @printf(io, "Data and simulated moments\n")
            @printf(io, "%-20s %16s %16s %16s\n", "Moment", "Data", "Model", "Residual")
            for jj in 1:length(momdiff)
                @printf(io, "%-20s %16.6f %16.6f %16.6f\n", 
                       momname[jj], datamoms[jj], simmoms[jj], momdiff[jj])
            end
            @printf(io, "\nParameter values\n")
            pnames = ["nu", "F_d", "kappa", "chi", "F_t"]
            for jj in eachindex(p)
                @printf(io, "%-10s %16.6f\n", pnames[jj], p[jj])
            end
            @printf(io, "\nGMM objective = %16.8f\n", bigQ)
            close(io)
        end
    else
        bigQ = 1e10  # Penalty for failed simulation
    end
    
    return bigQ::Float64
end


# ==========================================================================
# SMM statistics and inference
# ==========================================================================
function smmstats(p::Vector{Float64})
    pea = buildparam(p)
    
    simmoms = collect(momentgen(pea))
    datamoms = vec(collect(readdlm(kst.MOMS_FILE)))
    sw = collect(readdlm(kst.W_FILE))
    momname = vec(collect(readdlm(kst.MNAME_FILE)))
    pname = vec(collect(readdlm(kst.PNAME_FILE)))
    
    ch = sz.pick
    sw = sw[ch, ch]
    datamoms = datamoms[ch]
    simmoms = simmoms[ch]
    momname = momname[ch]
    
    if settings.complicated
        w = I(size(sw, 1))
    else
        w = inv(sw)
    end
    
    # Numerical gradient
    gee = grad(p, length(ch), sz.noestp)
    
    # Variance-covariance
    gwg = gee' * w * gee
    igwg = inv(gwg)
    
    wsw = w * sw * w
    gwswg = gee' * wsw * gee
    
    vc = igwg * gwswg * igwg
    
    # Adjustment for simulation error
    samplesize = 1000  # Approximate sample size
    nsim = float((sz.nYears - sz.burnin) * sz.nFirms) / float(samplesize)
    vc = (1.0 + 1.0/nsim) * vc
    
    standarderror = sqrt.(diag(vc))
    
    # J-test
    momdiff = datamoms .- simmoms
    jtest = (momdiff' * w * momdiff)[1]
    
    outtuple = (
        datamoms = datamoms,
        simmoms = simmoms,
        standarderror = standarderror,
        gee = gee,
        jtest = jtest,
        momdiff = momdiff,
        momname = momname,
        pname = pname
    )
    
    return outtuple
end


# ==========================================================================
# Numerical gradient
# ==========================================================================
function grad(x0, n, k)
    g = zeros(n, k)
    
    ax0 = abs.(x0)
    dax0 = sign.(x0)
    dh = 0.05 .* max.(ax0, 0.01) .* dax0
    
    for i in 1:k
        x_up = copy(x0)
        x_dw = copy(x0)
        x_up[i] += dh[i]
        x_dw[i] -= dh[i]
        
        g_up = fcn(x_up)
        g_dw = fcn(x_dw)
        
        g[:, i] = (g_up .- g_dw) ./ (2.0 * dh[i])
    end
    
    return g
end


# ==========================================================================
# Print results
# ==========================================================================
function print_smm_results(smm_results, p)
    datamoms = smm_results.datamoms
    simmoms = smm_results.simmoms
    standarderror = smm_results.standarderror
    gee = smm_results.gee
    jtest = smm_results.jtest
    momname = smm_results.momname
    pname = smm_results.pname
    
    io = open(kst.RESULTS_FILE, "w")
    
    @printf(io, "=" ^ 60 * "\n")
    @printf(io, "4D MODEL SMM RESULTS\n")
    @printf(io, "=" ^ 60 * "\n\n")
    
    @printf(io, "Parameter Estimates (Standard Errors)\n")
    @printf(io, "-" ^ 40 * "\n")
    for jj in eachindex(p)
        @printf(io, "%-10s = %10.4f (%8.4f)\n", pname[jj], p[jj], standarderror[jj])
    end
    
    @printf(io, "\n\nMoment Fit\n")
    @printf(io, "-" ^ 60 * "\n")
    @printf(io, "%-20s %12s %12s %12s\n", "Moment", "Data", "Model", "% Error")
    for jj in eachindex(datamoms)
        pct_err = 100 * abs(datamoms[jj] - simmoms[jj]) / max(abs(datamoms[jj]), 0.001)
        @printf(io, "%-20s %12.4f %12.4f %12.1f%%\n", 
               momname[jj], datamoms[jj], simmoms[jj], pct_err)
    end
    
    @printf(io, "\n\nJ-statistic = %.4f (df = %d)\n", jtest, length(datamoms) - length(p))
    
    @printf(io, "\n\nJacobian Matrix\n")
    @printf(io, "-" ^ 60 * "\n")
    @printf(io, "%-20s", "")
    for jj in 1:sz.noestp
        @printf(io, " %10s", pname[jj])
    end
    @printf(io, "\n")
    for ii in eachindex(momname)
        @printf(io, "%-20s", momname[ii])
        for jj in 1:sz.noestp
            @printf(io, " %10.4f", gee[ii, jj])
        end
        @printf(io, "\n")
    end
    
    close(io)
    
    println("\nResults saved to ", kst.RESULTS_FILE)
end
