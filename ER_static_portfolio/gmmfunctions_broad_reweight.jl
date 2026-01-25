# ==========================================================================
# 4D MODEL: Grid search / SMM estimation (FINAL, parallel-safe, prio-weighted)
#   - Uses prioritized W for estimation (W_pick6_prio125_diag.csv)
#   - Uses Σ for SEs (Sigma_pick6.csv)  [rename your file accordingly]
#   - Unique PROGRESS/EST/RESULTS/CKPT per run via RUNID
# ==========================================================================

ENV["GKSwstype"]          = "nul"
ENV["OPENBLAS_NUM_THREADS"]= "1"
ENV["MKL_NUM_THREADS"]     = "1"
ENV["OMP_NUM_THREADS"]     = "1"

using Random, Statistics, Printf, LinearAlgebra, Serialization
using DelimitedFiles
using Dates

# Include model files
include("durable_mod.jl")
include("collectfunctions.jl")
include("gmmfunctions_broad_reweight.jl")

using Main.sz, Main.kst, Main.settings, Main.dtp, Main.globals

# -------------------------------
# Run id + per-run output files
# -------------------------------
const RUNID = get(ENV, "RUNID",
    Dates.format(now(), "yyyymmdd_HHMMSS") * "_pid$(getpid())")

const PROGRESS_FILE = joinpath(kst.DATA_DIR, "progress_$(RUNID).txt")
const EST_FILE      = joinpath(kst.DATA_DIR, "est_$(RUNID).txt")
const RESULTS_FILE  = joinpath(kst.DATA_DIR, "results_$(RUNID).txt")
const CKPT_FILE     = joinpath("Output", "smm_ckpt_4D_prio125_$(RUNID).jls")

@inline time_left() = TDEAD - time()

# -------------------------------
# Files for estimation vs SEs
# -------------------------------
const WFILE_EST   = joinpath(kst.DATA_DIR, "W_pick6_prio125_diag.csv")  # weight used in objective
const SIGMA_FILE  = joinpath(kst.DATA_DIR, "Sigma_pick6.csv")           # covariance of moments for SEs (RECOMMENDED NAME)

# If you kept the old name, uncomment this instead:
# const SIGMA_FILE = joinpath(kst.DATA_DIR, "W_pick6_baseSigma.csv")

# ==========================================================================
# Load data moments and weighting matrix (estimation)
# ==========================================================================
const DAT = let
    dm     = vec(collect(readdlm(kst.MOMS_FILE)))
    W0     = collect(readdlm(WFILE_EST))
    mnames = vec(collect(readdlm(kst.MNAME_FILE)))
    pnames = vec(collect(readdlm(kst.PNAME_FILE)))
    (; dm, W0, mnames, pnames)
end

const CH   = sz.pick
const DM   = DAT.dm[CH]
const WRAW = DAT.W0[CH, CH]
const MN   = DAT.mnames[CH]
const PN   = DAT.pnames

const W = Symmetric((WRAW + WRAW')/2)   # already a weight matrix
const BIGPEN = 1e12

# ==========================================================================
# Estimation settings
# ==========================================================================
const BUDGET_MIN = 2860
const T0 = time()
const TDEAD = T0 + 60.0 * BUDGET_MIN
const seed = 1924
const n_trials = 2000

println("=" ^ 60)
println("4D MODEL ESTIMATION (RUNID=$(RUNID))")
println("=" ^ 60)
println("Host: $(gethostname())  Threads: $(Threads.nthreads())")
Random.seed!(seed)

# ==========================================================================
# Parameter bounds (4D model): [nu, F_d, kappa, chi, F_t]
# ==========================================================================
x_start = zeros(sz.noestp)
x_start[1] = 0.58
x_start[2] = 0.20
x_start[3] = 0.0001
x_start[4] = 0.50
x_start[5] = 0.10

lb = zeros(sz.noestp)
ub = zeros(sz.noestp)
lb[1] = 0.40;  ub[1] = 0.60
lb[2] = 0.001; ub[2] = 0.30
lb[3] = 0.0000;ub[3] = 0.005
lb[4] = 0.40;  ub[4] = 0.70
lb[5] = 0.001; ub[5] = 0.30

println("\nParameter bounds:")
pnames = ["nu", "F_d", "kappa", "chi", "F_t"]
for (i, name) in enumerate(pnames)
    @printf("  %s: [%.4f, %.4f], start = %.4f\n", name, lb[i], ub[i], x_start[i])
end

# ==========================================================================
# Checkpoint save/load (per-run)
# ==========================================================================
function save_ckpt!(x_best, f_best; path=CKPT_FILE)
    isdir(dirname(path)) || mkpath(dirname(path))
    serialize(path, (; x_best=copy(x_best), f_best=f_best, t=time()))
end

function load_ckpt(path=CKPT_FILE)
    isfile(path) ? deserialize(path) : nothing
end

# ==========================================================================
# Objective plumbing
# ==========================================================================

# --- parameter vector mapping ---
function buildparam(p::Vector{Float64})
    pea = ptrue(sz.nop)
    pea[5]  = p[1]   # nu
    pea[7]  = p[2]   # F_d
    pea[11] = p[3]   # kappa
    pea[16] = p[4]   # chi
    pea[17] = p[5]   # F_t
    return pea
end

# --- simulation wrapper (returns simulated moments, length == length(CH)) ---
function fcn(p::Vector{Float64})
    pea  = buildparam(p)
    moms = momentgen(pea)  # must return the picked moments in same order as DM
    bad = .!isfinite.(moms) .| (moms .== -100.0)
    if any(bad)
        @warn "Non-finite/flagged sim moments" first_bad=findfirst(bad) first_bad_val=moms[findfirst(bad)]
    end
    return moms
end

# --- GMM objective ---
function fcn(p::Vector{Float64}, fopt::Float64)
    simmoms = fcn(p)
    if !all(isfinite, simmoms)
        return BIGPEN
    end

    diff = DM .- simmoms
    bigQ = (diff' * W * diff)[1]
    if !isfinite(bigQ)
        return BIGPEN
    end

    # write progress only when improving
    if bigQ < fopt
        open(PROGRESS_FILE, "w") do io
            @printf(io, "RUNID: %s\n\n", RUNID)
            @printf(io, "Data vs. Simulated Moments (improved Q)\n\n")
            for j in eachindex(diff)
                @printf(io, "%-20s  data = %12.6f   sim = %12.6f   resid = %12.6f\n",
                        MN[j], DM[j], simmoms[j], diff[j])
            end
            @printf(io, "\nParameters:\n")
            for j in eachindex(p)
                @printf(io, "%-20s  %12.6f\n", PN[j], p[j])
            end
            @printf(io, "\nGMM objective Q = %12.6f\n", bigQ)
        end
        open(EST_FILE, "w") do io
            @printf(io, "RUNID: %s\n", RUNID)
            for j in eachindex(p)
                @printf(io, "%25.16f\n", p[j])
            end
        end
    end
    return bigQ
end

# --- safe wrapper for search loop ---
function safe_fcn(x, best_so_far)
    try
        return fcn(x, best_so_far)
    catch e
        @warn "fcn crashed" x exception=(e, catch_backtrace())
        return Inf
    end
end

# ==========================================================================
# Gradient + SEs (uses Σ and the same W used in estimation)
# ==========================================================================

function grad(x0::Vector{Float64}, n::Int, k::Int)
    g = zeros(n, k)
    ax0  = abs.(x0)
    dax0 = sign.(x0) .+ (x0 .== 0.0)
    dh   = 1e-3 .* max.(ax0, 1e-2) .* dax0

    x = copy(x0)
    for i in 1:k
        xi = x0[i]

        x[i] = xi + dh[i]
        mup  = fcn(x)
        x[i] = xi - dh[i]
        mdw  = fcn(x)
        x[i] = xi

        if !all(isfinite, mup) || !all(isfinite, mdw)
            g[:, i] .= 0.0
        else
            g[:, i] = (mup .- mdw) ./ (2.0 * dh[i])
        end
    end
    return g
end

@inline function _safe_pinv(A; rtol=1e-8, atol=0.0, ridge=0.0)
    A = Symmetric((A + A')/2)
    if ridge > 0
        return inv(A + ridge*I)
    end
    S = svd(Matrix(A); full=false)
    tol = max(atol, rtol * maximum(S.S))
    r = count(>(tol), S.S)
    if r == 0
        return zeros(size(A))
    end
    Sinv = Diagonal(vcat(1 ./ S.S[1:r], zeros(length(S.S)-r)))
    return S.Vt' * Sinv * S.U'
end

function smmstats(p::Vector{Float64}; n_sample::Int=10_000)
    simmoms = fcn(p)
    if !all(isfinite, simmoms)
        error("smmstats: simulated moments are not finite; cannot compute SEs.")
    end

    datamoms = vec(collect(readdlm(kst.MOMS_FILE)))[CH]
    Σfull    = collect(readdlm(SIGMA_FILE))
    Σ        = Symmetric((Σfull[CH, CH] + Σfull[CH, CH]')/2)

    # weight matrix used in estimation (already loaded as W)
    W_est = W

    diff = datamoms .- simmoms
    n = length(CH); k = length(p)

    G = grad(p, n, k)

    GWG  = G' * W_est * G
    igwg = _safe_pinv(GWG; ridge=1e-10)

    vc = igwg * (G' * W_est * Σ * W_est * G) * igwg
    vc = Symmetric((vc + vc')/2)
    vc = Matrix(vc) .* (1.0 + 1.0/n_sample)

    se = sqrt.(abs.(diag(vc)))
    jtest = (diff' * W_est * diff)[1]

    return (datamoms=datamoms, simmoms=simmoms, se=se, jacobian=G, vcov=vc, jtest=jtest)
end

function print_smm_results(stats, p::Vector{Float64})
    open(RESULTS_FILE, "w") do io
        @printf(io, "RUNID: %s\n\n", RUNID)
        @printf(io, "Final SMM results (nu, F_d, kappa, chi, F_t)\n\n")
        @printf(io, "Parameter estimates and standard errors\n\n")
        for j in eachindex(p)
            @printf(io, "%-20s  %12.6f   (SE = %12.6f)\n", PN[j], p[j], stats.se[j])
        end

        @printf(io, "\nData vs. Simulated Moments and residuals\n\n")
        for j in eachindex(stats.datamoms)
            resid = stats.datamoms[j] - stats.simmoms[j]
            @printf(io, "%-20s  data = %12.6f   sim = %12.6f   resid = %12.6f\n",
                    MN[j], stats.datamoms[j], stats.simmoms[j], resid)
        end

        @printf(io, "\nOver-identification (J-test): %12.6f\n", stats.jtest)
        @printf(io, "(df = %d)\n", length(stats.datamoms) - length(p))
    end
end

# ==========================================================================
# Run: initial point + checkpoint + random search
# ==========================================================================
println("\nTesting initial point...")
f_test = safe_fcn(x_start, BIGPEN)
println("Initial objective: ", f_test)

ck = load_ckpt()
if ck === nothing || length(ck.x_best) != sz.noestp
    x_opt = copy(x_start)
    f_opt = safe_fcn(x_opt, BIGPEN)
    save_ckpt!(x_opt, f_opt)
else
    x_opt = copy(ck.x_best)
    f_opt = ck.f_best
    println("Resumed from checkpoint: f = ", f_opt)
end

println("\nStarting search...")
println("x0 = ", x_opt)
println("f(x0) = ", f_opt)

global nevals = 0
t_loop0 = time()

for i in 1:n_trials
    if time_left() <= 30.0
        break
    end

    x = lb .+ (ub .- lb) .* rand(sz.noestp)
    f = safe_fcn(x, f_opt)
    nevals += 1

    if isfinite(f) && f < f_opt
        x_opt = copy(x)
        f_opt = f
        @printf("\n*** NEW OPTIMUM at trial %d ***\n", i)
        @printf("f = %.6g\n", f_opt)
        for (j, name) in enumerate(pnames)
            @printf("  %s = %.6f\n", name, x_opt[j])
        end
        save_ckpt!(x_opt, f_opt)
    end

    if (i % 25 == 0) || (time_left() < 60.0)
        save_ckpt!(x_opt, f_opt)
        elapsed = time() - t_loop0
        @printf("Trial %d/%d, evals=%d, %.2fs/eval, time_left=%.0fs, best_f=%.2e\n",
                i, n_trials, nevals, elapsed/max(nevals,1), time_left(), f_opt)
    end
end

println("\n" * "=" ^ 60)
println("Search complete (RUNID=$(RUNID))")
println("=" ^ 60)
println("Evaluations: ", nevals)
println("Best objective: ", f_opt)
println("Best parameters:")
for (j, name) in enumerate(pnames)
    @printf("  %s = %.6f\n", name, x_opt[j])
end

save_ckpt!(x_opt, f_opt)

if isfinite(f_opt)
    println("\nComputing SMM standard errors...")
    stats = smmstats(x_opt)
    print_smm_results(stats, x_opt)
else
    @warn "Best objective is not finite; skipping standard errors."
end

println("\nWrote:")
println("  ", PROGRESS_FILE)
println("  ", EST_FILE)
println("  ", RESULTS_FILE)
println("  ", CKPT_FILE)
