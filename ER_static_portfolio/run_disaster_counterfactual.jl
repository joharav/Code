#!/usr/bin/env julia
using Random, DelimitedFiles, Distributions, LinearAlgebra, Plots, Statistics, Printf
using StatsBase, KernelDensity, JLD2, PrettyTables, DataFrames, CSV

include("durable_mod.jl")
include("collectfunctions.jl")
include("disaster_process.jl")   # provides makegrids_disaster
using Main.sz, Main.settings

p = ptrue(sz.nop)

# baseline builder
gridA = makegrids

# scenarios
pi_annual = 0.2
kappas    = [log(1.18), log(1.89), log(2.5)]
labels    = ["mild", "moderate", "severe"]

println("=== Baseline (no disaster) ===")
m_base = momentgen(p; grid_builder=gridA)
println("Baseline moments:\n", m_base)

nm = length(m_base)
ns = length(kappas)
store = zeros(ns, 4 + nm)

for (ik, κ) in enumerate(kappas)
    label = labels[ik]
    println("\n=== Disaster scenario: $label ===")
    println("π^y = $(pi_annual), κ_e_log = $(κ)  (≈ $(exp(κ))x devaluation)")

    gridB = (q::Vector{Float64}) -> begin
        g_dis, _ = makegrids_disaster(q; pi_annual=pi_annual, kappa_e_log=κ)
        g_dis
    end

    # welfare (A vs B) without touching globals
    res_w = welfare_summary(p, p; gridA=gridA, gridB=gridB)
    println("Across steady states (A→B): ", res_w.acrossSS, " %")
    println("CEV_BA (B vs A): ", 100*res_w.cev_BA, " %")

    # moments in B
    m_dis = momentgen(p; grid_builder=gridB)
    println("Disaster moments:\n", m_dis)

    store[ik, 1] = pi_annual
    store[ik, 2] = κ
    store[ik, 3] = res_w.acrossSS
    store[ik, 4] = res_w.cev_BA
    store[ik, 5:end] .= m_dis
end

header = ["pi_annual","kappa_e_log","acrossSS","cev_BA",
          "duration_mean","dwealth_mean","dwealth_var","adj_rate","usd_share_mean","usd_share_vol",
          "cons_vol","d_spend_vol","a_eff_vol"]

outpath = "Output/disaster_counterfactuals.txt"
isdir(dirname(outpath)) || mkpath(dirname(outpath))

open(outpath, "w") do io
    println(io, join(header, '\t'))
    writedlm(io, store)
end

println("\nSaved disaster summary table to $(outpath)")
