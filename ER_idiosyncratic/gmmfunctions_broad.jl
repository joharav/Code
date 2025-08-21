using LinearAlgebra, DelimitedFiles, Printf, Statistics
using Main.kst

# -------------------------------
# === Parameter vector utilities ===
# -------------------------------

function buildparam(p::Vector{Float64})
    pea = ptrue(sz.nop)
    pea[5]  = p[1]   # ν_ndurables
    pea[7]  = p[2]   # f_d
    pea[11] = p[3]   # f_t
    return pea
end

# -------------------------------
# === Simulation wrapper ===
# -------------------------------

function fcn(p::Vector{Float64})
    pea  = buildparam(p)
    moms = momentgen(pea)
    return moms[sz.pick]   # pick only d_income_ratio, d_wealth_ratio, adj_ratio
end

# -------------------------------
# === GMM objective ===
# -------------------------------
function fcn(p::Vector{Float64}, fopt::Float64)
    simmoms = fcn(p)
    if simmoms[1] == -100.0
        return sz.maxfunc
    end

    datamoms = vec(collect(readdlm(kst.MOMS_FILE)))
    Wraw     = collect(readdlm(kst.W_FILE))
    momname  = collect(readdlm(kst.MNAME_FILE))
    pname    = collect(readdlm(kst.PNAME_FILE))

    # Subset to selected moments
    ch       = sz.pick
    datamoms = datamoms[ch]
    Wraw     = Wraw[ch, ch]
    momname  = momname[ch]

    # Weighting matrix
    W = settings.complicated ? I(size(Wraw,1)) : inv(Wraw)

    # Objective
    diff = datamoms .- simmoms
    bigQ = (diff' * W * diff)[1]

    # Progress report if improved
    if bigQ < fopt
        open(kst.PROGRESS_FILE, "w") do io
            @printf(io, "Data vs. Simulated Moments (improved Q)\n\n")
            for j in eachindex(diff)
                @printf(io, "%-20s  data = %12.6f   sim = %12.6f\n",
                        momname[j], datamoms[j], simmoms[j])
            end
            @printf(io, "\nParameters:\n")
            for j in eachindex(p)
                @printf(io, "%-20s  %12.6f\n", pname[j], p[j])
            end
            @printf(io, "\nGMM objective Q = %12.6f\n", bigQ)
        end

        # Save last best parameters
        open(kst.EST_FILE, "w") do io
            for j in eachindex(p); @printf(io, "%25.16f\n", p[j]); end
        end
    end

    return bigQ
end

# -------------------------------
# === Numerical gradient ===
# -------------------------------
function grad(x0::Vector{Float64}, n::Int, k::Int)
    g = zeros(n, k)
    ax0  = abs.(x0)
    dax0 = sign.(x0) .+ (x0 .== 0.0)
    dh   = 1e-3 .* max.(ax0, 1e-2) .* dax0

    xdup = x0 .+ dh
    xddw = x0 .- dh

    argup = repeat(x0', k, 1)
    argdw = repeat(x0', k, 1)

    for i in 1:k
        argup[i,i] = xdup[i]
        argdw[i,i] = xddw[i]
    end

    g_up = zeros(n,k)
    g_dw = zeros(n,k)
    for i in 1:k
        g_up[:,i] = fcn(view(argup,i,:))
        g_dw[:,i] = fcn(view(argdw,i,:))
    end

    for i in 1:k
        g[:,i] = (g_up[:,i] .- g_dw[:,i]) ./ (2.0 * dh[i])
    end
    return g
end

# -------------------------------
# === SMM statistics ===
# -------------------------------
function smmstats(p::Vector{Float64}; n_sample::Int=10_000)
    simmoms = fcn(p)
    datamoms = vec(collect(readdlm(kst.MOMS_FILE)))
    Σ        = collect(readdlm(kst.W_FILE))
    momname  = collect(readdlm(kst.MNAME_FILE))

    ch       = sz.pick
    datamoms = datamoms[ch]
    Σ        = Σ[ch, ch]
    momname  = momname[ch]

    W = settings.complicated ? I(size(Σ,1)) : inv(Σ)

    diff = datamoms .- simmoms
    n = length(ch)
    k = length(p)
    G = grad(p, n, k)

    GWG  = G' * W * G
    igwg = inv(GWG)
    vc   = igwg * (G' * W * Σ * W * G) * igwg
    vc  .*= (1.0 + 1.0/n_sample)

    se = sqrt.(diag(vc))
    jtest = (diff' * W * diff)[1]

    return (datamoms=datamoms, simmoms=simmoms, se=se,
            jacobian=G, vcov=vc, jtest=jtest)
end

# -------------------------------
# === Save SMM results ===
# -------------------------------
function print_smm_results(stats, p::Vector{Float64})
    momname = collect(readdlm(kst.MNAME_FILE))[sz.pick]
    pname   = collect(readdlm(kst.PNAME_FILE))

    open(kst.RESULTS_FILE, "w") do io
        @printf(io, "Final SMM results (ν_ndurables, f_d, f_t)\n\n")
        @printf(io, "Parameter estimates and standard errors\n\n")
        for j in eachindex(p)
            @printf(io, "%-20s  %12.6f   (SE = %12.6f)\n", pname[j], p[j], stats.se[j])
        end

        @printf(io, "\nData vs. Simulated Moments and residuals\n\n")
        for j in eachindex(stats.datamoms)
            resid = stats.datamoms[j] - stats.simmoms[j]
            @printf(io, "%-28s  data = %12.6f   sim = %12.6f   resid = %12.6f\n",
                    momname[j], stats.datamoms[j], stats.simmoms[j], resid)
        end

        @printf(io, "\nOver-identification (J-test): %12.6f\n", stats.jtest)
        @printf(io, "(df = %d)\n", length(stats.datamoms) - length(p))
    end
end
