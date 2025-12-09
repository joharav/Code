# gridsearch_durables.jl — robust, time-budgeted random search

ENV["GKSwstype"] = "nul"
ENV["OPENBLAS_NUM_THREADS"] = "1"
ENV["MKL_NUM_THREADS"]     = "1"
ENV["OMP_NUM_THREADS"]      = "1"


@inline time_left() = TDEAD - time()

using Random, Statistics, Printf, LinearAlgebra, Serialization, Distributions
# using JLD2, DataFrames, CSV, PrettyTables, StatsBase, KernelDensity, Plots

include("durable_mod.jl")
include("collectfunctions.jl")
include("simann.jl")
include("gmmfunctions_broad.jl")

using Main.sz, Main.kst, Main.settings, Main.dtp, Main.globals


# -------------------------------
# === Data & weights (cached) ===
# -------------------------------
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
const ALL_MN = vec(collect(readdlm(kst.MNAME_FILE)))  # full names (pre-pick)
const BUDGET_MIN = 2860                    # in-Julia budget; leaves cushion
const T0 = time()
const TDEAD = T0 + 60.0 * BUDGET_MIN
const EPS  = 1e-12
const seed = 1924
const n_trials = 2000

# ---------- robust wrappers ----------
function safe_fcn(x, best_so_far)
    try
        return fcn(x, best_so_far)   # fcn calls buildparam internally
    catch e
        @warn "fcn crashed" x exception=(e, catch_backtrace())
        return Inf
    end
end

function save_ckpt!(x_best, f_best; path="Output/smm_ckpt.jls")
    isdir(dirname(path)) || mkpath(dirname(path))
    serialize(path, (; x_best=copy(x_best), f_best=f_best, t=time()))
end

function load_ckpt(path="Output/smm_ckpt.jls")
    isfile(path) ? deserialize(path) : nothing
end

# ----------------- Script body -----------------
println("Host: $(gethostname())  Threads: $(Threads.nthreads())")
Random.seed!(seed)

# Bounds & start
x_start = zeros(sz.noestp)
x_start[1] = 0.473848   # nu
x_start[2] = 0.050269   # f_d
x_start[3] = 0.824595   # kappa
x_start[4] = 0.497794     # chi
x_start[5] = 0.320361       # ft

lb = zeros(sz.noestp);  ub = zeros(sz.noestp)
lb[1] = 0.4;  ub[1] = 0.6
lb[2] = 0.01; ub[2] = 0.1
lb[3] = 0.7; ub[3] = 2.0
lb[4] = 0.35;  ub[4] = 0.7
lb[5] = 0.2;  ub[5] = 0.8

try
    _ = safe_fcn(x_start, 1e12)
catch
end

# Resume from checkpoint if present
ck = load_ckpt()
if ck === nothing || length(ck.x_best) != sz.noestp
    @warn "Ignoring incompatible or missing checkpoint" ck_len = (ck === nothing ? missing : length(ck.x_best)) noestp = sz.noestp
    x_opt = copy(x_start)
    f_opt = safe_fcn(x_opt, 1e12)
    save_ckpt!(x_opt, f_opt)  
else
    x_opt = copy(ck.x_best)
    f_opt = ck.f_best
end
println("Starting search with x0 = ", x_opt, ", f(x0) = ", f_opt)


global i = 1
global nevals = 0
t_loop0 = time()
while i ≤ n_trials && time_left() > 30.0      
    x = lb .+ (ub .- lb) .* rand(sz.noestp)
    f = safe_fcn(x, f_opt)
    global nevals += 1
    if isfinite(f) && f < f_opt
        global x_opt = copy(x)   
        global f_opt = f        
        @printf("New optimum at trial %d: f = %.6g, x = %s\n", i, f_opt, string(x_opt))
        save_ckpt!(x_opt, f_opt)
    end
    # save occasionally regardless, and near the end
    if (i % 25 == 0) || (time_left() < 60.0)
        save_ckpt!(x_opt, f_opt)
        elapsed = time() - t_loop0
        @printf("Progress: i=%d, evals=%d, avg %.3fs/eval, time_left=%.1fs\n",
                i, nevals, elapsed / max(nevals,1), time_left())
    end
    global i += 1
end

println("\n⏱️ Stopped at i=$i with time_left=$(round(time_left(); digits=1))s")
println("Best so far: f = ", f_opt, " at x = ", x_opt)
save_ckpt!(x_opt, f_opt)

if isfinite(f_opt)
    stats = smmstats(x_opt)
    print_smm_results(stats, x_opt)
else
    @warn "Best objective is not finite; skipping SMM standard errors."
end
