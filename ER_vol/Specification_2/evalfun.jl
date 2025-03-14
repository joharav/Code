using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2;
include("../durable_mod.jl");
include("collectfunctions.jl");

using Main.sz, Main.kst, Main.settings, Main.globals, Main.dtp; 

pea_2 = ptrue(sz.nop); 

#answ=valfun(pea);   

moms_2 = momentgen(pea);
#include("compstat.jl");


println("done with solving the model")
#include("compare_utilities.jl");