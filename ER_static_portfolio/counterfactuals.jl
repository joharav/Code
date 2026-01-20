# ==========================================================================
# 4D MODEL: Counterfactual analysis (robust)
# ==========================================================================

using DelimitedFiles, Statistics, Printf, LinearAlgebra
using Random, Distributions

include("durable_mod.jl")
include("collectfunctions.jl")
include("ptrue.jl")
include("welfare.jl")

using Main.sz, Main.kst, Main.settings, Main.globals

# ----------------------------
# Canonical parameter indices
# ----------------------------
const P = (;
    beta=1, delta=2, rho_e=3, sigma_e=4, nu=5, gamma=6,
    F_d=7, wage=8, r_foreign=9, pd=10, kappa=11, tau=12, h=13,
    rho_y=14, sigma_y=15,
    chi=16, ft=17
)

# ----------------------------
# Est load + build baseline
# ----------------------------
function load_theta_hat()
    est_file = joinpath(kst.OUT_DIR, "estfil.txt")
    if isfile(est_file)
        θ = vec(collect(readdlm(est_file)))
        return Float64.(θ)
    else
        # MUST match buildparam(θ) internal expectations
        return [0.55, 1.5, 0.10, 0.50, 0.30]  # example: nu, F_d, kappa, chi, ft
    end
end

function baseline_params()
    θhat = load_theta_hat()
    p = buildparam(θhat)

    # sanity checks
    @assert length(p) >= P.ft "pea length smaller than expected index map"
    @assert p[P.beta] > 0 && p[P.beta] < 1
    @assert p[P.sigma_e] >= 0
    @assert p[P.kappa] >= 0
    return p
end

# ----------------------------
# Counterfactual modifiers
# ----------------------------
cf_higher_Fd(p; factor=2.0) = (q=copy(p); q[P.F_d] *= factor; q)
cf_higher_sigma_e(p; new_sigma=0.30) = (q=copy(p); q[P.sigma_e] = new_sigma; q)
cf_fixed_ER(p) = (q=copy(p); q[P.sigma_e] = 0.0; q)
cf_higher_delta(p; new_delta=0.05) = (q=copy(p); q[P.delta] = new_delta; q)
cf_lower_kappa(p; new_kappa=0.01) = (q=copy(p); q[P.kappa] = new_kappa; q)

# This is *not* a hard "no dollar saving" unless you also restrict s choices.
cf_high_kappa(p; kappa=10.0) = (q=copy(p); q[P.kappa] = kappa; q)

# ----------------------------
# Moments unpacker (guarded)
# ----------------------------
function unpack_moments(m::AbstractVector)
    @assert length(m) >= 6 "momentgen output changed; update unpack_moments."
    return (;
        duration_mean = m[1],
        dwealth_mean  = m[2],
        dwealth_var   = m[3],
        adj_rate      = m[4],
        dollar_share  = m[5],
        dollar_vol    = m[6],
    )
end

# ----------------------------
# Welfare wrapper (robust)
# ----------------------------
function welfare_change(pe_A::Vector{Float64}, pe_B::Vector{Float64})
    res = welfare_full_summary(pe_A, pe_B)
    return res.cev_BA
end

# ----------------------------
# OPTIONAL: hard "no dollar saving"
# Only works if your solver reads grids from parameters/settings.
# The clean way is: rebuild grids with s_grid = [0.0] and re-solve.
# ----------------------------
function momentgen_nodollar_hard(p::Vector{Float64})
    # If your code supports it, temporarily override the s-grid.
    # This assumes you have something like `build_grids(p)` or `solve_model(p; grids=...)`.
    # If you DON'T, delete this function and stick to high-kappa as a soft proxy.

    g = build_grids(p)                 # <-- you must have this or equivalent
    g = merge(g, (s = [0.0],))         # force only s=0 choice
    answ = solve_model(p; grids=g)     # <-- you must have this or equivalent
    sim  = simmodel(answ, p)           # or whatever momentgen uses internally
    return makemoments(sim, p)         # return same vector format as momentgen
end

# ----------------------------
# Main runner (deterministic)
# ----------------------------
function run_counterfactuals()
    p_base = baseline_params()

    println("="^60)
    println("4D MODEL COUNTERFACTUAL ANALYSIS")
    println("="^60)

    # Baseline
    println("\n=== BASELINE ===")
    moms_base = momentgen(p_base)
    mech_base = unpack_moments(moms_base)

    @printf("duration_mean  = %.4f years\n", mech_base.duration_mean)
    @printf("adj_rate       = %.4f\n", mech_base.adj_rate)
    @printf("dwealth_mean   = %.4f\n", mech_base.dwealth_mean)
    @printf("dollar_share   = %.4f\n", mech_base.dollar_share)

    scenarios = [
        ("higher_Fd",        p -> cf_higher_Fd(p; factor=2.0),         "Higher F^d (×2)"),
        ("fixed_ER",         p -> cf_fixed_ER(p),                       "Fixed ER (σ_e=0)"),
        ("higher_sigma_e",   p -> cf_higher_sigma_e(p; new_sigma=0.30), "Higher σ_e (0.30)"),
        ("lower_kappa",      p -> cf_lower_kappa(p; new_kappa=0.01),    "Lower κ (0.01)"),
        ("high_kappa_soft",  p -> cf_high_kappa(p; kappa=10.0),         "High κ (10) [soft no-dollar]"),
        # ("no_dollar_hard", p -> p,                                     "No dollar [hard s=0]")  # only if momentgen_nodollar_hard works
    ]

    results = Dict{String, NamedTuple}()

    for (key, modfun, label) in scenarios
        println("\n--- CF: $label ---")
        p_cf = modfun(p_base)

        moms_cf = momentgen(p_cf)
        mech_cf = unpack_moments(moms_cf)

        λ = try
            welfare_change(p_base, p_cf)
        catch
            NaN
        end

        @printf("Δduration_mean = %+.4f\n", mech_cf.duration_mean - mech_base.duration_mean)
        @printf("Δadj_rate      = %+.4f\n", mech_cf.adj_rate - mech_base.adj_rate)
        @printf("Δdollar_share  = %+.4f\n", mech_cf.dollar_share - mech_base.dollar_share)
        @printf("λ (welfare)    = %.4f%%\n", λ)

        results[key] = (moms=mech_cf, welfare=λ)
    end

    # Summary table (stable order)
    println("\n" * "="^60)
    println("SUMMARY")
    println("="^60)

    @printf("\n%-20s %10s %10s %10s %10s\n", "Scenario", "Adj Rate", "Duration", "\$ Share", "Welfare")
    @printf("%-20s %10.4f %10.4f %10.4f %10s\n",
            "Baseline", mech_base.adj_rate, mech_base.duration_mean, mech_base.dollar_share, "-")

    for (key, _, _) in scenarios
        res = results[key]
        @printf("%-20s %10.4f %10.4f %10.4f %+10.2f%%\n",
                key, res.moms.adj_rate, res.moms.duration_mean, res.moms.dollar_share, res.welfare)
    end

    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_counterfactuals()
end
