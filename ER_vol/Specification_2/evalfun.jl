using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2;
include("durable_mod.jl");
include("collectfunctions.jl");

using Main.sz, Main.settings, Main.globals, Main.dtp; 

pea = ptrue(sz.nop); 

#answ=valfun(pea);   

moms_2, answ_fixed_er = momentgen(pea);
#include("compstat.jl");


println("done with solving the model")
#include("compare_utilities.jl");