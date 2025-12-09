# counterfactuals.jl

using DelimitedFiles, Statistics, Printf, LinearAlgebra
using Random, Distributions, Plots, StatsBase, KernelDensity, JLD2, PrettyTables, DataFrames, CSV;

include("durable_mod.jl")
include("collectfunctions.jl")
include("momentgen.jl")
include("makemoments_5assets.jl")    
include("ptrue.jl")
include("welfare_cases.jl")   

using Main.sz, Main.kst, Main.settings, Main.globals


function load_theta_hat()
    est_file = joinpath(kst.OUT_DIR, "estfil.txt")
    θ = vec(collect(readdlm(est_file)))
    return Float64.(θ)
end


function baseline_params()
    θhat = load_theta_hat()
    return buildparam(θhat)
end

# ---- counterfactual parameter tweaks ----

# indices (from your ptrue):
# 1 beta, 2 delta, 3 rho_e, 4 sigma_e, 5 nu, 6 gamma, 7 F^d,
# 8 w, 9 r_foreign, 10 pd, 11 kappa_a, 12 tau, 13 h,
# 14 rho_y, 15 sigma_y, 16 chi, 17 ft

function cf_higher_Fd(p::Vector{Float64}; factor=2.0)
    q = copy(p)
    q[7] *= factor
    return q
end

function cf_higher_sigma_e(p::Vector{Float64}; factor=2.0)
    q = copy(p)
    q[4] *= factor
    return q
end

function cf_fixed_ER(p::Vector{Float64})
    q = copy(p)
    q[4] = 0.0    # sigma_e = 0
    return q
end

function cf_no_dollar_saving(p::Vector{Float64})
    q = copy(p)
    q[11] = 100.0   # kappa_a huge
    q[17] = 100.0   # ft huge if that enters foreign asset choice
    return q
end

# ---- helper to compute welfare change A→B ----

function welfare_change(pe_A::Vector{Float64}, pe_B::Vector{Float64})
    res = welfare_full_summary(pe_A, pe_B)
    return res.λ_composite
end
"""
    unpack_moments(m::AbstractVector)

Map moment vector to a NamedTuple with the semantic names used below.
Ordering must match makemoments_5assets.jl.
"""
function unpack_moments(m::AbstractVector)
    @assert length(m) >= 9 "Expected at least 9 moments, got $(length(m))"
    return (
        duration_mean   = m[1],
        dwealth_mean    = m[2],
        dwealth_var     = m[3],
        adj_rate        = m[4],
        usd_share_mean  = m[5],
        usd_share_var   = m[6],
        cons_vol        = m[7],
        d_spend_vol     = m[8],
        a_eff_vol       = m[9],
    )
end

# ---- main driver ----

function run_counterfactuals()
    p_base = baseline_params()

    println("=== Baseline ===")
    mech_base = unpack_moments(momentgen(p_base))
   
    @printf("duration_mean       = %.4f\n", mech_base.duration_mean)
    @printf("adj_rate            = %.4f\n", mech_base.adj_rate)
    @printf("owner_share         = %.4f\n", mech_base.dwealth_mean)   # if you dropped owner_share, adjust text
    @printf("usd_share_mean      = %.4f\n", mech_base.usd_share_mean)
    @printf("cons_vol            = %.4f\n", mech_base.cons_vol)
    @printf("d_spend_vol         = %.4f\n", mech_base.d_spend_vol)
    @printf("a_eff_vol           = %.4f\n", mech_base.a_eff_vol)

    # ---- higher F^d ----
    println("\n--- CF: higher F^d ---")
    p_Fd     = cf_higher_Fd(p_base; factor=2.0)
    mech_Fd  = unpack_moments(momentgen(p_Fd))
    λ_Fd     = welfare_change(p_base, p_Fd)

    @printf("Δduration_mean      = %.4f\n", mech_Fd.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate           = %.4f\n", mech_Fd.adj_rate        - mech_base.adj_rate)
    @printf("Δdwealth_var        = %.4f\n", mech_Fd.dwealth_var     - mech_base.dwealth_var)
    @printf("Δusd_share_mean     = %.4f\n", mech_Fd.usd_share_mean  - mech_base.usd_share_mean)
    @printf("Δcons_vol           = %.4f\n", mech_Fd.cons_vol        - mech_base.cons_vol)
    @printf("Δd_spend_vol        = %.4f\n", mech_Fd.d_spend_vol     - mech_base.d_spend_vol)
    @printf("Δa_eff_vol          = %.4f\n", mech_Fd.a_eff_vol       - mech_base.a_eff_vol)
    @printf("λ (higher F^d)      = %.4f\n", λ_Fd)

    # ---- higher sigma_e ----
    println("\n--- CF: higher σ_e ---")
    p_se    = cf_higher_sigma_e(p_base; factor=2.0)
    mech_se = unpack_moments(momentgen(p_se))
    λ_se    = welfare_change(p_base, p_se)

    @printf("Δduration_mean      = %.4f\n", mech_se.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate           = %.4f\n", mech_se.adj_rate        - mech_base.adj_rate)
    @printf("Δdwealth_var        = %.4f\n", mech_se.dwealth_var     - mech_base.dwealth_var)
    @printf("Δusd_share_mean     = %.4f\n", mech_se.usd_share_mean  - mech_base.usd_share_mean)
    @printf("Δcons_vol           = %.4f\n", mech_se.cons_vol        - mech_base.cons_vol)
    @printf("Δd_spend_vol        = %.4f\n", mech_se.d_spend_vol     - mech_base.d_spend_vol)
    @printf("Δa_eff_vol          = %.4f\n", mech_se.a_eff_vol       - mech_base.a_eff_vol)
    @printf("λ (higher σ_e)      = %.4f\n", λ_se)

    # ---- fixed ER ----
    println("\n--- CF: fixed ER ---")
    p_fix    = cf_fixed_ER(p_base)
    mech_fix = unpack_moments(momentgen(p_fix))
    λ_fix    = welfare_change(p_base, p_fix)

    @printf("Δduration_mean      = %.4f\n", mech_fix.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate           = %.4f\n", mech_fix.adj_rate        - mech_base.adj_rate)
    @printf("Δdwealth_var        = %.4f\n", mech_fix.dwealth_var     - mech_base.dwealth_var)
    @printf("Δusd_share_mean     = %.4f\n", mech_fix.usd_share_mean  - mech_base.usd_share_mean)
    @printf("Δcons_vol           = %.4f\n", mech_fix.cons_vol        - mech_base.cons_vol)
    @printf("Δd_spend_vol        = %.4f\n", mech_fix.d_spend_vol     - mech_base.d_spend_vol)
    @printf("Δa_eff_vol          = %.4f\n", mech_fix.a_eff_vol       - mech_base.a_eff_vol)
    @printf("λ (fixed ER)        = %.4f\n", λ_fix)

    # ---- no dollar saving ----
    println("\n--- CF: no dollar saving ---")
    p_nusd    = cf_no_dollar_saving(p_base)
    mech_nusd = unpack_moments(momentgen(p_nusd))
    λ_nusd    = welfare_change(p_base, p_nusd)

    @printf("Δduration_mean      = %.4f\n", mech_nusd.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate           = %.4f\n", mech_nusd.adj_rate        - mech_base.adj_rate)
    @printf("Δdwealth_var        = %.4f\n", mech_nusd.dwealth_var     - mech_base.dwealth_var)
    @printf("Δusd_share_mean     = %.4f\n", mech_nusd.usd_share_mean  - mech_base.usd_share_mean)
    @printf("Δcons_vol           = %.4f\n", mech_nusd.cons_vol        - mech_base.cons_vol)
    @printf("Δd_spend_vol        = %.4f\n", mech_nusd.d_spend_vol     - mech_base.d_spend_vol)
    @printf("Δa_eff_vol          = %.4f\n", mech_nusd.a_eff_vol       - mech_base.a_eff_vol)
    @printf("λ (no dollar)       = %.4f\n", λ_nusd)

    return nothing
end

run_counterfactuals()

