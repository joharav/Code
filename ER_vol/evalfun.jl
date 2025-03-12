using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2;
include("durable_mod.jl");
include("collectfunctions.jl");

using Main.sz, Main.kst, Main.settings, Main.globals, Main.dtp; 

pea = ptrue(sz.nop); 

#answ=valfun(pea);   
#simdata = simmodel(answ);
#simdata_irf = simmodel_girf(answ, Int(sz.nYears/2));
#girf = girf_plots(simdata_irf, simdata);
#cirf_c = compute_cirf(vec(girf[1]), 8, "c");
#cirf_d = compute_cirf(vec(girf[2]), 8, "d");
#cirf_a = compute_cirf(vec(girf[3]), 8, "a");;
moms = momentgen(pea);
#include("compstat.jl");
