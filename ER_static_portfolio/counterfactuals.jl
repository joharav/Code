# ==========================================================================
# 4D MODEL: Counterfactual analysis
# ==========================================================================

using DelimitedFiles, Statistics, Printf, LinearAlgebra
using Random, Distributions

include("durable_mod.jl")
include("collectfunctions.jl")
include("ptrue.jl")
include("welfare.jl")

using Main.sz, Main.kst, Main.settings, Main.globals


function load_theta_hat()
    est_file = joinpath(kst.OUT_DIR, "estfil.txt")
    if isfile(est_file)
        θ = vec(collect(readdlm(est_file)))
        return Float64.(θ)
    else
        # Use default starting values
        return [0.55, 1.5, 0.10, 0.50, 0.30]  # nu, F_d, kappa, chi, F_t
    end
end

function baseline_params()
    θhat = load_theta_hat()
    return buildparam(θhat)
end

# ==========================================================================
# Counterfactual parameter modifications
# Indices: 1=beta, 2=delta, 3=rho_e, 4=sigma_e, 5=nu, 6=gamma, 7=F_d,
#          8=wage, 9=r_foreign, 10=pd, 11=kappa, 12=tau, 13=h,
#          14=rho_y, 15=sigma_y, 16=chi, 17=F_t
# ==========================================================================

# Higher durable fixed cost
function cf_higher_Fd(p::Vector{Float64}; factor=2.0)
    q = copy(p)
    q[7] *= factor  # F_d
    return q
end

# Higher exchange rate volatility
function cf_higher_sigma_e(p::Vector{Float64}; new_sigma=0.30)
    q = copy(p)
    q[4] = new_sigma
    return q
end

# Fixed exchange rate (no volatility)
function cf_fixed_ER(p::Vector{Float64})
    q = copy(p)
    q[4] = 0.0  # sigma_e = 0
    return q
end

# No dollar saving (very high transaction cost)
function cf_no_dollar_saving(p::Vector{Float64})
    q = copy(p)
    q[11] = 10.0  # kappa very high
    return q
end

# Lower dollar transaction cost
function cf_lower_kappa(p::Vector{Float64}; new_kappa=0.01)
    q = copy(p)
    q[11] = new_kappa
    return q
end

# Higher depreciation
function cf_higher_delta(p::Vector{Float64}; new_delta=0.05)
    q = copy(p)
    q[2] = new_delta
    return q
end


# ==========================================================================
# Unpack moments helper
# ==========================================================================
function unpack_moments(m::AbstractVector)
    @assert length(m) >= 6 "Expected at least 6 moments"
    return (
        duration_mean = m[1],
        dwealth_mean = m[2],
        dwealth_var = m[3],
        adj_rate = m[4],
        dollar_share = m[5],
        dollar_vol = m[6],
    )
end


# ==========================================================================
# Welfare change computation
# ==========================================================================
function welfare_change(pe_A::Vector{Float64}, pe_B::Vector{Float64})
    res = welfare_full_summary(pe_A, pe_B)
    return res.cev_BA
end


# ==========================================================================
# Main counterfactual analysis
# ==========================================================================
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
    
    results = Dict{String, NamedTuple}()
    
    # ---------------------------------------------------------------------
    # CF1: Higher F^d (doubles durable adjustment cost)
    # ---------------------------------------------------------------------
    println("\n--- CF1: Higher F^d (×2) ---")
    p_Fd = cf_higher_Fd(p_base; factor=2.0)
    moms_Fd = momentgen(p_Fd)
    mech_Fd = unpack_moments(moms_Fd)
    
    try
        λ_Fd = welfare_change(p_base, p_Fd)
    catch
        λ_Fd = NaN
    end
    
    @printf("Δduration_mean = %+.4f\n", mech_Fd.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate      = %+.4f\n", mech_Fd.adj_rate - mech_base.adj_rate)
    @printf("Δdollar_share  = %+.4f\n", mech_Fd.dollar_share - mech_base.dollar_share)
    @printf("λ (welfare)    = %.4f%%\n", λ_Fd)
    
    results["higher_Fd"] = (moms=mech_Fd, welfare=λ_Fd)
    
    # ---------------------------------------------------------------------
    # CF2: Fixed exchange rate
    # ---------------------------------------------------------------------
    println("\n--- CF2: Fixed ER (σ_e = 0) ---")
    p_fix = cf_fixed_ER(p_base)
    moms_fix = momentgen(p_fix)
    mech_fix = unpack_moments(moms_fix)
    
    try
        λ_fix = welfare_change(p_base, p_fix)
    catch
        λ_fix = NaN
    end
    
    @printf("Δduration_mean = %+.4f\n", mech_fix.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate      = %+.4f\n", mech_fix.adj_rate - mech_base.adj_rate)
    @printf("Δdollar_share  = %+.4f\n", mech_fix.dollar_share - mech_base.dollar_share)
    @printf("λ (welfare)    = %.4f%%\n", λ_fix)
    
    results["fixed_ER"] = (moms=mech_fix, welfare=λ_fix)
    
    # ---------------------------------------------------------------------
    # CF3: No dollar saving
    # ---------------------------------------------------------------------
    println("\n--- CF3: No dollar saving (κ = 10) ---")
    p_nusd = cf_no_dollar_saving(p_base)
    moms_nusd = momentgen(p_nusd)
    mech_nusd = unpack_moments(moms_nusd)
    
    try
        λ_nusd = welfare_change(p_base, p_nusd)
    catch
        λ_nusd = NaN
    end
    
    @printf("Δduration_mean = %+.4f\n", mech_nusd.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate      = %+.4f\n", mech_nusd.adj_rate - mech_base.adj_rate)
    @printf("Δdollar_share  = %+.4f\n", mech_nusd.dollar_share - mech_base.dollar_share)
    @printf("λ (welfare)    = %.4f%%\n", λ_nusd)
    
    results["no_dollar"] = (moms=mech_nusd, welfare=λ_nusd)
    
    # ---------------------------------------------------------------------
    # CF4: Higher ER volatility
    # ---------------------------------------------------------------------
    println("\n--- CF4: Higher σ_e (0.30) ---")
    p_se = cf_higher_sigma_e(p_base; new_sigma=0.30)
    moms_se = momentgen(p_se)
    mech_se = unpack_moments(moms_se)
    
    try
        λ_se = welfare_change(p_base, p_se)
    catch
        λ_se = NaN
    end
    
    @printf("Δduration_mean = %+.4f\n", mech_se.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate      = %+.4f\n", mech_se.adj_rate - mech_base.adj_rate)
    @printf("Δdollar_share  = %+.4f\n", mech_se.dollar_share - mech_base.dollar_share)
    @printf("λ (welfare)    = %.4f%%\n", λ_se)
    
    results["higher_sigma_e"] = (moms=mech_se, welfare=λ_se)
    
    # ---------------------------------------------------------------------
    # CF5: Lower kappa (easier dollar access)
    # ---------------------------------------------------------------------
    println("\n--- CF5: Lower κ (0.01) ---")
    p_lk = cf_lower_kappa(p_base; new_kappa=0.01)
    moms_lk = momentgen(p_lk)
    mech_lk = unpack_moments(moms_lk)
    
    try
        λ_lk = welfare_change(p_base, p_lk)
    catch
        λ_lk = NaN
    end
    
    @printf("Δduration_mean = %+.4f\n", mech_lk.duration_mean - mech_base.duration_mean)
    @printf("Δadj_rate      = %+.4f\n", mech_lk.adj_rate - mech_base.adj_rate)
    @printf("Δdollar_share  = %+.4f\n", mech_lk.dollar_share - mech_base.dollar_share)
    @printf("λ (welfare)    = %.4f%%\n", λ_lk)
    
    results["lower_kappa"] = (moms=mech_lk, welfare=λ_lk)
    
    # ---------------------------------------------------------------------
    # Summary table
    # ---------------------------------------------------------------------
    println("\n" * "="^60)
    println("SUMMARY")
    println("="^60)
    
    @printf("\n%-20s %10s %10s %10s %10s\n", 
           "Scenario", "Adj Rate", "Duration", "$ Share", "Welfare")
    @printf("%-20s %10.4f %10.2f %10.4f %10s\n",
           "Baseline", mech_base.adj_rate, mech_base.duration_mean, 
           mech_base.dollar_share, "-")
    
    for (name, res) in results
        @printf("%-20s %10.4f %10.2f %10.4f %+10.2f%%\n",
               name, res.moms.adj_rate, res.moms.duration_mean,
               res.moms.dollar_share, res.welfare)
    end
    
    return results
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_counterfactuals()
end
