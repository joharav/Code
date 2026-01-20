#!/usr/bin/env julia
using Random, DelimitedFiles, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2, PrettyTables, DataFrames, CSV

include("durable_mod.jl")   # defines grid_builder = p -> makegrids(p)
include("collectfunctions.jl")
using Main.sz, Main.settings

# ------------------------------------------------------------------
# Baseline parameter vector
# ------------------------------------------------------------------
p = ptrue(sz.nop)

println("=== Baseline (no disaster) ===")
m_base = momentgen(p)
println("Baseline moments (duration_mean, dwealth_mean, dwealth_var, adj_rate, usd_share_mean, usd_share_vol, cons_vol, d_spend_vol, a_eff_vol):")
println(m_base)

# ------------------------------------------------------------------
# Crisis scenarios: κ_e ∈ {ln 1.2, ln 1.4, ln 1.6}, π^y ≈ 2%
# ------------------------------------------------------------------
pi_annual = 0.1                    # annual disaster probability
kappas    = [log(1.2), log(1.5), log(2.0)]
labels    = ["mild", "moderate", "severe"]

nm = length(m_base)
ns = length(kappas)

# store [pi_annual, kappa_e, acrossSS, cev_BA, m1..m9] per scenario
store = zeros(ns, 4 + nm)

for (ik, κ) in enumerate(kappas)
    label = labels[ik]
    println("\n=== Disaster scenario: $label ===")
    println("π^y = $(pi_annual), κ_e = $(κ)  (≈ $(exp(κ))x devaluation)")

    # Welfare comparison (baseline grids vs disaster grids, same p)
    res_w = welfare_disaster(p; pi_annual = pi_annual, kappa_e = κ)

    println("Across steady states (A→B, baseline→disaster): ", res_w.acrossSS, " %")
    println("CEV_BA (B vs A, composite): ", 100 * res_w.cev_BA, " %")

    # Moments under disaster process:
    # Temporarily point the grid builder to disaster grids
    local grid_builder = q -> makegrids_disaster(q; pi_annual = pi_annual, kappa_e = κ)
    m_dis = momentgen(p)
    println("Disaster moments (duration_mean, dwealth_mean, dwealth_var, adj_rate, usd_share_mean, usd_share_vol, cons_vol, d_spend_vol, a_eff_vol):")
    println(m_dis)

    # store row
    store[ik, 1] = pi_annual
    store[ik, 2] = κ
    store[ik, 3] = res_w.acrossSS
    store[ik, 4] = res_w.cev_BA
    store[ik, 5:end] .= m_dis

end

# ------------------------------------------------------------------
# Dump scenarios table to disk
# ------------------------------------------------------------------
header = ["pi_annual",
          "kappa_e",
          "acrossSS",
          "cev_BA",
          "duration_mean",
          "dwealth_mean",
          "dwealth_var",
          "adj_rate",
          "usd_share_mean",
          "usd_share_vol",
          "cons_vol",
          "d_spend_vol",
          "a_eff_vol"]

outpath = "Output/disaster_counterfactuals.txt"
isdir(dirname(outpath)) || mkpath(dirname(outpath))

open(outpath, "w") do io
    println(io, join(header, '\t'))   # or ',' if you want CSV
    writedlm(io, store)              # numeric matrix is fine
end


println("\nSaved disaster summary table to $(outpath)")
