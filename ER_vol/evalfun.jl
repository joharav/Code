using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2;
include("durable_mod.jl");
include("collectfunctions.jl");

using Main.sz,  Main.settings, Main.globals, Main.dtp; 

pea = ptrue(sz.nop); 

moms, answ_baseline = momentgen(pea);

if settings.specif_three
    include("Specification_3/evalfun.jl")

end

if settings.specif_two
    include("Specification_2/evalfun.jl")
end

# Store policy objects for comparison
policies = Dict(
    "baseline" => answ_baseline,
    "high_vol" => answ_high_vol,
    "fixed_er" => answ_fixed_er
)
# Plot policy functions for comparisonfixed_er

plot_policy_functions(policies)

if settings.compstat
    include("compstat.jl")
end

println("done with solving the model")


