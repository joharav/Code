# Function to compute Consumption Equivalent Variation (CEV)
function compute_cev(v_pre::Vector{Float64}, v_post::Vector{Float64}, ppp::Vector{Float64})
    nu      = ppp[5]     # Nondurable share parameter
    gamma   = ppp[6]  # Risk-aversion coefficient

    ratio   = mean(v_post) / mean(v_pre)
    cev     = (ratio)^(1 / ((1 - gamma) * (1 - nu))) - 1
    return cev * 100  # Convert to percentage
end

# Function to compute dispersion measures
function compute_dispersion(x::Vector{Float64})
    std_dev = std(x)
    iqr     = percentile(x, 75) - percentile(x, 25)
    p90_10  = percentile(x, 90) - percentile(x, 10)
    
    outuple = (std_dev, iqr, p90_10)
    return outuple
end
function compute_cev_array(v_alt::Array{Float64,4}, v_base::Vector{Float64}, pe::Vector{Float64})
    γ = pe[6]
    if γ == 1.0
        cev = exp.(v_alt .- reshape(v_base, 1, 1, 1, :)) .- 1
    else
        cev = ((v_alt ./ reshape(v_base, 1, 1, 1, :)) .^ (1 / (1 - γ))) .- 1
    end
    return cev
end
