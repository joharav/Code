# ==========================================================================
# 4D MODEL: Welfare computations
# ==========================================================================

using Statistics
using StatsBase: percentile

# Compute Consumption Equivalent Variation (CEV)
function compute_cev(v_pre::Vector{Float64}, v_post::Vector{Float64}, ppp::Vector{Float64})
    nu = ppp[5]       # Non-durable share
    gamma = ppp[6]    # Risk aversion
    
    ratio = mean(v_post) / mean(v_pre)
    cev = (ratio)^(1 / ((1 - gamma) * (1 - nu))) - 1
    return cev * 100  # percentage
end

# CEV from value arrays
function compute_cev(v_pre::Array{Float64,4}, v_post::Array{Float64,4}, ppp::Vector{Float64})
    nu = ppp[5]
    gamma = ppp[6]
    
    ratio = mean(v_post) / mean(v_pre)
    cev = (ratio)^(1 / ((1 - gamma) * (1 - nu))) - 1
    return cev * 100
end

# CEV array: compare alternative to baseline at each state
function compute_cev_array(v_alt::Array{Float64,4}, v_base::Array{Float64,4}, pe::Vector{Float64})
    γ = pe[6]
    
    if γ == 1.0
        cev = exp.(v_alt .- v_base) .- 1
    else
        cev = ((v_alt ./ max.(v_base, 1e-10)) .^ (1 / (1 - γ))) .- 1
    end
    return cev * 100  # percentage
end

# Dispersion measures
function compute_dispersion(x::Vector{Float64})
    x_fin = x[isfinite.(x)]
    std_dev = std(x_fin)
    iqr = percentile(x_fin, 75) - percentile(x_fin, 25)
    p90_10 = percentile(x_fin, 90) - percentile(x_fin, 10)
    
    return (std_dev=std_dev, iqr=iqr, p90_10=p90_10)
end

# Full welfare summary between baseline and counterfactual
function welfare_full_summary(pe_base::Vector{Float64}, pe_alt::Vector{Float64})
    # Solve both economies
    ans_base = valfun(pe_base)
    ans_alt = valfun(pe_alt)
    
    # CEV from value functions
    cev_BA = compute_cev(ans_base.v, ans_alt.v, pe_base)
    
    # Ergodic distributions
    dist_base = compute_ergodic(ans_base)
    dist_alt = compute_ergodic(ans_alt)
    
    # Expected value under distributions
    EV_base = sum(ans_base.v .* dist_base)
    EV_alt_keepdist = sum(ans_alt.v .* dist_base)  # Alt value, base distribution
    EV_alt = sum(ans_alt.v .* dist_alt)             # Alt value, alt distribution
    
    # Welfare changes
    keepDistAB = (EV_alt_keepdist / EV_base)^(1/((1-pe_base[6])*(1-pe_base[5]))) - 1
    acrossSS = (EV_alt / EV_base)^(1/((1-pe_base[6])*(1-pe_base[5]))) - 1
    
    return (
        cev_BA = cev_BA,
        λ_composite = cev_BA,
        keepDistAB = keepDistAB * 100,
        acrossSS = acrossSS * 100,
        EV_base = EV_base,
        EV_alt = EV_alt
    )
end
