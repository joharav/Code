using Random, Distributions, Statistics, Printf, DelimitedFiles, LinearAlgebra, Plots, JLD2, 
    StatsBase, KernelDensity, DataFrames, CSV, PrettyTables;

# ----------------- Include your functions -----------------
include("durable_mod.jl")
include("collectfunctions.jl")
include("simann.jl")
include("gmmfunctions_broad.jl")

using Main.sz
using Main.kst
using Main.settings
using Main.dtp
using Main.globals


# ----------------- Flags -----------------
restart = false          # restart from previous estimation
doiterations = true      # run optimization or just compute stats

# ----------------- Starting values and bounds -----------------
x_start = zeros(sz.noestp)
x_start[1] = 0.6   # nu
x_start[2] = 0.4   # f_d
x_start[3] = 0.4   # f_t

lb = zeros(sz.noestp)
ub = zeros(sz.noestp)

lb[1] = 0.3;   ub[1] = 0.999   # nu
lb[2] = 0.01;   ub[2] = 0.9     # f_d
lb[3] = 0.01;   ub[3] = 0.9     # f_t
pea = buildparam(x_start)
moms = momentgen(pea)
# ----------------- SA tuning -----------------
iseed = 1924
n_trials = 2000
Random.seed!(iseed)
# ---------------------------
# Simulated Annealing / Random Search
# ---------------------------

# Initialize best guess
x_opt = copy(x_start)                  # start with user-defined x_start
f_opt = fcn(x_opt, 1e12)              # evaluate f at x_start
println("Starting search with x_start = ", x_start, ", f(x_start) = ", f_opt)


# Perform random search for remaining trials
for i in 2:n_trials
    # Random draw within bounds for each parameter
    x = lb .+ (ub .- lb) .* rand(sz.noestp)

    f = fcn(x, f_opt)                  # evaluate objective

    if f < f_opt
        global x_opt = copy(x)
        global f_opt = f
        println("New optimum found at trial $i: x = ", x_opt, ", f = ", f_opt)
    end
end

println("\n✅ Optimization finished")
println("Optimal parameters: ", x_opt)
println("Function value at optimum: ", f_opt)
