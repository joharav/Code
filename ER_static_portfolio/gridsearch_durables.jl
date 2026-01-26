# ==========================================================================
# 4D MODEL: Grid search / SMM estimation
# ==========================================================================

ENV["GKSwstype"] = "nul"
ENV["OPENBLAS_NUM_THREADS"] = "1"
ENV["MKL_NUM_THREADS"]     = "1"
ENV["OMP_NUM_THREADS"]      = "1"

@inline time_left() = TDEAD - time()

using Random, Statistics, Printf, LinearAlgebra, Serialization, Distributions
using DelimitedFiles

# Include model files
include("durable_mod.jl")
include("collectfunctions.jl")
include("gmmfunctions_broad.jl")

using Main.sz, Main.kst, Main.settings, Main.dtp, Main.globals


# ==========================================================================
# Load data moments and weighting matrix
# ==========================================================================
const DAT = let
    dm = vec(collect(readdlm(kst.MOMS_FILE)))
    W0 = collect(readdlm(kst.W_FILE))
    mnames = vec(collect(readdlm(kst.MNAME_FILE)))
    pnames = vec(collect(readdlm(kst.PNAME_FILE)))
    (; dm, W0, mnames, pnames)
end

const CH   = sz.pick
const DM   = DAT.dm[CH]
const WRAW = DAT.W0[CH, CH]
const MN   = DAT.mnames[CH]
const PN   = DAT.pnames
const W    = let Wsym = Symmetric((WRAW + WRAW')/2); ridge=1e-10; inv(Wsym + ridge*I); end
const BIGPEN = 1e12

# Time budget
const BUDGET_MIN = 2860                    
const T0 = time()
const TDEAD = T0 + 60.0 * BUDGET_MIN
const seed = 1924
const n_trials = 2000


# ==========================================================================
# Safe wrapper for objective function
# ==========================================================================
function safe_fcn(x, best_so_far)
    try
        return fcn(x, best_so_far)
    catch e
        @warn "fcn crashed" x exception=(e, catch_backtrace())
        return Inf
    end
end


# ==========================================================================
# Checkpoint save/load
# ==========================================================================
function save_ckpt!(x_best, f_best; path="Output/smm_ckpt_4D.jls")
    isdir(dirname(path)) || mkpath(dirname(path))
    serialize(path, (; x_best=copy(x_best), f_best=f_best, t=time()))
end

function load_ckpt(path="Output/smm_ckpt_4D.jls")
    isfile(path) ? deserialize(path) : nothing
end


# ==========================================================================
# Main estimation script
# ==========================================================================
println("=" ^ 60)
println("4D MODEL ESTIMATION")
println("=" ^ 60)
println("Host: $(gethostname())  Threads: $(Threads.nthreads())")
Random.seed!(seed)

# Parameter bounds - UPDATED for 4D model
# Order: [nu, F_d, kappa, chi, F_t]
x_start = zeros(sz.noestp)
x_start[1] = 0.5       # nu (non-durable share)
x_start[2] = 0.2        # F_d (durable fixed cost) - INCREASED
x_start[3] = 0.0001       # kappa (dollar transaction cost) - DECREASED

lb = zeros(sz.noestp)
ub = zeros(sz.noestp)

# Updated bounds based on moment mismatch diagnosis:
# - F_d needs to be MUCH higher to reduce adjustment rate
# - kappa needs to be LOWER to increase dollar share
lb[1] = 0.20;  ub[1] = 0.60   # nu
lb[2] = 0.001;  ub[2] = 0.7   # F_d - WIDER, HIGHER range
lb[3] = 0.0000;  ub[3] = 0.009   # kappa - LOWER range

println("\nParameter bounds:")
pnames = ["nu", "F_d", "kappa"]
for (i, name) in enumerate(pnames)
    @printf("  %s: [%.2f, %.2f], start = %.2f\n", name, lb[i], ub[i], x_start[i])
end

# Test initial point
println("\nTesting initial point...")
try
    f_test = safe_fcn(x_start, 1e12)
    println("Initial objective: ", f_test)
catch e
    println("Initial test failed: ", e)
end

# Resume from checkpoint if present
ck = load_ckpt()
if ck === nothing || length(ck.x_best) != sz.noestp
    if ck !== nothing
        @warn "Ignoring incompatible checkpoint"
    end
    x_opt = copy(x_start)
    f_opt = safe_fcn(x_opt, 1e12)
    save_ckpt!(x_opt, f_opt)
else
    x_opt = copy(ck.x_best)
    f_opt = ck.f_best
    println("Resumed from checkpoint: f = ", f_opt)
end

println("\nStarting search...")
println("x0 = ", x_opt)
println("f(x0) = ", f_opt)


# ==========================================================================
# Random search with local refinement
# ==========================================================================
global i = 1
global nevals = 0
t_loop0 = time()

while i ≤ n_trials && time_left() > 30.0
    # Random point in bounds
    x = lb .+ (ub .- lb) .* rand(sz.noestp)
    
    f = safe_fcn(x, f_opt)
    global nevals += 1
    
    if isfinite(f) && f < f_opt
        global x_opt = copy(x)
        global f_opt = f
        @printf("\n*** NEW OPTIMUM at trial %d ***\n", i)
        @printf("f = %.6g\n", f_opt)
        for (j, name) in enumerate(pnames)
            @printf("  %s = %.4f\n", name, x_opt[j])
        end
        save_ckpt!(x_opt, f_opt)
    end
    
    # Progress report
    if (i % 25 == 0) || (time_left() < 60.0)
        save_ckpt!(x_opt, f_opt)
        elapsed = time() - t_loop0
        @printf("Trial %d/%d, evals=%d, %.2fs/eval, time_left=%.0fs, best_f=%.2e\n",
                i, n_trials, nevals, elapsed/max(nevals,1), time_left(), f_opt)
    end
    
    global i += 1
end

println("\n" * "=" ^ 60)
println("Search complete")
println("=" ^ 60)
println("Trials: ", i-1)
println("Evaluations: ", nevals)
println("Best objective: ", f_opt)
println("Best parameters:")
for (j, name) in enumerate(pnames)
    @printf("  %s = %.6f\n", name, x_opt[j])
end

save_ckpt!(x_opt, f_opt)

# Compute standard errors if converged
if isfinite(f_opt)
    println("\nComputing SMM standard errors...")
    try
        stats = smmstats(x_opt)
        print_smm_results(stats, x_opt)
    catch e
        @warn "Standard error computation failed" exception=(e, catch_backtrace())
    end
else
    @warn "Best objective is not finite; skipping standard errors."
end
