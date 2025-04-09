using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2, DataStructures;
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
policies = OrderedDict(
    "Baseline" => answ_baseline,
    "High Volatility" => answ_high_vol,
    "Fixed Exchange Rate" => answ_fixed_er
)
# Plot policy functions for comparisonfixed_er

plot_policy_functions(policies)

if settings.compstat
    include("compstat.jl")
end

println("done with solving the model")


