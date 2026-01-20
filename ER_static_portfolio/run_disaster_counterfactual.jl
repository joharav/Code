#!/usr/bin/env julia
using Random, DelimitedFiles, Distributions, LinearAlgebra, Plots, Statistics, Printf
using StatsBase, KernelDensity, JLD2, PrettyTables, DataFrames, CSV

include("durable_mod.jl")         # should define: makegrids, maybe grid_builder, momentgen uses it
include("collectfunctions.jl")    # momentgen, welfare_disaster, etc.
using Main.sz, Main.settings

# ----------------------------
# Baseline parameters
# ----------------------------
p = ptrue(sz.nop)

# ----------------------------
# Helper: disaster grid builder (returns ONLY grids)
# ----------------------------
function grid_builder_disaster(q::Vector{Float64}; pi_annual::Float64, kappa_e_log::Float64)
    g_dis, _meta = makegrids_disaster(q; pi_annual=pi_annual, kappa_e_log=kappa_e_log)
    return g_dis
end

# ----------------------------
# Ensure baseline grid builder is set (momentgen must use Main.grid_builder)
# ----------------------------
Main.grid_builder = makegrids   # baseline

println("=== Baseline (no disaster) ===")
m_base = momentgen(p)
println("Baseline moments:")
println(m_base)

# ----------------------------
# Crisis scenarios
# ----------------------------
pi_annual = 0.10
kappas    = [log(1.2), log(1.5), log(2.0)]
labels    = ["mild", "moderate", "severe"]

nm = length(m_base)
ns = length(kappas)

# store: [pi_annual, kappa_e_log, acrossSS, cev_BA, moments...]
store = zeros(ns, 4 + nm)

for (ik, κ) in enumerate(kappas)
    label = labels[ik]
    println("\n=== Disaster scenario: $label ===")
    println("π^y = $(pi_annual), κ_e = $(κ)  (≈ $(exp(κ))x devaluation)")

    # 1) Welfare comparison (must internally solve both cases with correct grid builders)
    #    This call is only valid if welfare_disaster itself uses makegrids_disaster OR accepts a grid_builder.
    #    If welfare_disaster uses Main.grid_builder, then set it before calling welfare_disaster too.
    Main.grid_builder = q -> grid_builder_disaster(q; pi_annual=pi_annual, kappa_e_log=κ)

    res_w = welfare_disaster(p; pi_annual=pi_annual, kappa_e_log=κ)

    println("Across steady states (A→B, baseline→disaster): ", res_w.acrossSS, " %")
    println("CEV_BA (B vs A, composite): ", 100 * res_w.cev_BA, " %")

    # 2) Moments under disaster process
    m_dis = momentgen(p)
    println("Disaster moments:")
    println(m_dis)

    # store row
    store[ik, 1] = pi_annual
    store[ik, 2] = κ
    store[ik, 3] = res_w.acrossSS
    store[ik, 4] = res_w.cev_BA
    store[ik, 5:end] .= m_dis
end

# Restore baseline builder
Main.grid_builder = makegrids

# ----------------------------
# Dump scenarios table
# ----------------------------
header = ["pi_annual",
          "kappa_e_log",
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
    println(io, join(header, '\t'))
    writedlm(io, store)
end

println("\nSaved disaster summary table to $(outpath)")
