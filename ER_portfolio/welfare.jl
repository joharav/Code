# Function to compute Consumption Equivalent Variation (CEV)
function compute_cev(v_pre::Vector{Float64}, v_post::Vector{Float64}, ppp::Vector{Float64})
    nu      = ppp[5]     # Nondurable share parameter
    gamma   = ppp[6]  # Risk-aversion coefficient

    ratio   = mean(v_post) / mean(v_pre)
    cev     = (ratio)^(1 / ((1 - gamma) * (1 - nu))) - 1
   return cev * 100  # Convert to percentage
end

# === 5D version for two-asset state v[ie,iy,iaa,ia,id] ===
function compute_cev_array(v_alt::Array{Float64,5}, v_base::Vector{Float64}, pe::Vector{Float64})
    γ = pe[6]
    # v_base is assumed along the last dim (id) like before
    if γ == 1.0
        cev = exp.(v_alt .- reshape(v_base, 1, 1, 1, 1, :)) .- 1
    else
        cev = ((v_alt ./ reshape(v_base, 1, 1, 1, 1, :)) .^ (1 / (1 - γ))) .- 1
    end
    return cev
end
