using LinearAlgebra, DelimitedFiles, Printf, Statistics
using Main.kst

# -------------------------------
# === Parameter vector utilities ===
# -------------------------------

function buildparam(p::Vector{Float64})
    pea = ptrue(sz.nop)
    pea[5]  = p[1]   # ν_ndurables
    pea[7]  = p[2]   # f_d, fixed cost
   # pea[11] = p[3]   # kappa, fixed cost asset
   # pea[16] = p[4]   # chi, maintenance effectiveness 
    #pea[17] = p[5]   # ft, time fixed cost
    return pea
end

# -------------------------------
# === Simulation wrapper (robust) ===
# -------------------------------


function fcn(p::Vector{Float64})
    pea = buildparam(p)
    moms= momentgen(pea)   
    bad = .!isfinite.(moms) .| (moms .== -100.0)
    if any(bad)
        @warn "Non-finite/flagged sim moments" first_bad=findfirst(bad) first_bad_val=moms[findfirst(bad)]
    end
    return moms  # if moms is already the 6 you need, just `return moms`
end



# -------------------------------
# === GMM objective (robust) ===
# -------------------------------
function fcn(p::Vector{Float64}, fopt::Float64)
    simmoms = fcn(p)                           # uses robust sim wrapper above
    if !all(isfinite, simmoms)
        return BIGPEN
    end

    diff = DM .- simmoms
    bigQ = (diff' * W * diff)[1]
    if !isfinite(bigQ)
        return BIGPEN
    end

    if bigQ < fopt
        open(kst.PROGRESS_FILE, "w") do io
            @printf(io, "Data vs. Simulated Moments (improved Q)\n\n")
            for j in eachindex(diff)
                @printf(io, "%-20s  data = %12.6f   sim = %12.6f\n", MN[j], DM[j], simmoms[j])
            end
            @printf(io, "\nParameters:\n")
            for j in eachindex(p)
                @printf(io, "%-20s  %12.6f\n", PN[j], p[j])
            end
            @printf(io, "\nGMM objective Q = %12.6f\n", bigQ)
        end
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
    dax0 = sign.(x0) .+ (x0 .== 0.0)   # avoids 0 step
    dh   = 1e-3 .* max.(ax0, 1e-2) .* dax0

    # Central differences without huge allocations
    x = copy(x0)
    for i in 1:k
        xi = x0[i]

        x[i] = xi + dh[i]
        mup  = fcn(x)
        x[i] = xi - dh[i]
        mdw  = fcn(x)
        x[i] = xi

        if !all(isfinite, mup) || !all(isfinite, mdw)
            g[:,i] .= 0.0   # if sim fails at perturbed points, neutral slope
        else
            g[:,i] = (mup .- mdw) ./ (2.0 * dh[i])
        end
    end
    return g
end


# tiny helper
@inline function _safe_pinv(A; rtol=1e-8, atol=0.0, ridge=0.0)
    m, n = size(A)
    if m==0 || n==0
        throw(ArgumentError("safe_pinv: empty matrix ($m×$n)"))
    end
    A = Symmetric((A + A')/2)  # for square GWG
    if ridge > 0
        return inv(A + ridge*I)
    end
    # fall back to SVD-based pseudo-inverse
    S = svd(Matrix(A); full=false)
    tol = max(atol, rtol * maximum(S.S))
    r = count(>(tol), S.S)
    if r == 0
        return zeros(n, m)
    end
    Sinv = Diagonal(vcat(1 ./ S.S[1:r], zeros(length(S.S)-r)))
    return S.Vt' * Sinv * S.U'
end

function smmstats(p::Vector{Float64}; n_sample::Int=10_000)
    simmoms = fcn(p)
    if !all(isfinite, simmoms)
        error("smmstats: simulated moments are not finite; cannot compute SEs.")
    end

    datamoms = vec(collect(readdlm(kst.MOMS_FILE)))
    Σ        = collect(readdlm(kst.W_FILE))
    momname  = collect(readdlm(kst.MNAME_FILE))

    ch       = sz.pick
    datamoms = datamoms[ch]
    Σ        = Σ[ch, ch]
    momname  = momname[ch]

    Σ = Symmetric((Σ + Σ')/2)
    Wridge = 1e-10
    W = settings.complicated ? I(size(Σ,1)) : inv(Σ + Wridge*I)

    diff = datamoms .- simmoms
    n = length(ch)
    k = length(p)
    @assert n >= 1 "No selected moments (n==0); check sz.pick"
    @assert k >= 1 "No parameters (k==0)"

    G = grad(p, n, k)
    if any(!isfinite, G)
        # If gradient fails anywhere, avoid SVD explosions
        se = fill(NaN, k)
        jtest = (diff' * W * diff)[1]
        return (datamoms=datamoms, simmoms=simmoms, se=se,
                jacobian=G, vcov=fill(NaN, k, k), jtest=jtest)
    end

    GWG = G' * W * G
    # ridge for stability of the inverse
    igwg = _safe_pinv(GWG; ridge=1e-10)

    vc   = igwg * (G' * W * Σ * W * G) * igwg
    vc   = (vc + vc')/2                 # symmetrize
    vc  .*= (1.0 + 1.0/n_sample)

    se = sqrt.(abs.(diag(vc)))          # abs for tiny negative diag from FP
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
        @printf(io, "Final SMM results (ν_ndurables, f_d, kappa)\n\n")
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
