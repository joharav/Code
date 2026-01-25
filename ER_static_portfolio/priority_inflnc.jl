using Statistics
using StatsBase: Weights, mean, var
using LinearAlgebra
using DataFrames, CSV
using DelimitedFiles: writedlm

include("durable_mod.jl")
using Main.kst

# ------------------------------------------------------------
# Helpers: finite/missing mask + aligned x,w (keeps weight scale)
# ------------------------------------------------------------
function finmask(xu::AbstractVector)
    m = falses(length(xu))
    @inbounds for i in eachindex(xu)
        if !ismissing(xu[i])
            xi = Float64(xu[i])
            m[i] = isfinite(xi)
        end
    end
    return m
end

function aligned_xw(xu::AbstractVector, w::Vector{Float64})
    mask = finmask(xu)
    x  = Float64.(xu[mask])
    ww = w[mask]                  # DO NOT renormalize here
    return x, ww, mask
end

# ------------------------------------------------------------
# Influence functions (return full-length vector with zeros off mask)
# ------------------------------------------------------------
function IF_wmean(xu::AbstractVector, w::Vector{Float64})
    x, ww, mask = aligned_xw(xu, w)
    μ = mean(x, Weights(ww))
    out = zeros(Float64, length(xu))
    out[mask] = x .- μ
    return out
end

function IF_wvar(xu::AbstractVector, w::Vector{Float64})
    x, ww, mask = aligned_xw(xu, w)
    μ  = mean(x, Weights(ww))
    σ2 = sum(ww .* (x .- μ).^2) / sum(ww)   # prob-weight variance
    out = zeros(Float64, length(xu))
    out[mask] = (x .- μ).^2 .- σ2
    return out
end

function IF_wshare_pos(xu::AbstractVector, w::Vector{Float64})
    x, ww, mask = aligned_xw(xu, w)
    ind = x .> 0.0
    θ = sum(ww .* ind) / sum(ww)
    out = zeros(Float64, length(xu))
    out[mask] = Float64.(ind) .- θ
    return out
end

# ------------------------------------------------------------
# Load microdata and data-moment vector (10-moment file)
# ------------------------------------------------------------
df = CSV.read(joinpath(kst.DATA_DIR, "EFHU_moments_data_weighted.csv"), DataFrame)

mrow_full = CSV.read(joinpath(kst.DATA_DIR, "moments_vector_data_weighted.csv"), DataFrame)
datamom_full = vec(Matrix(mrow_full)[1, :])   # length 10

# Pick 6 moments in this order (must match IF columns below)
const pick6 = [1, 3, 4, 6, 8, 5]              # duration, dwealth_mean, dwealth_var, adj_rate, usd_share_mean, usd_share_var
datamom = datamom_full[pick6]

mom_names = [
    "duration_mean",
    "dwealth_mean",
    "dwealth_var",
    "adj_rate",
    "dollar_share",
    "dollar_vol",
]

# ------------------------------------------------------------
# Weights
# ------------------------------------------------------------
N = nrow(df)
w_raw = Float64.(df.pesoEFHU)
w = w_raw ./ sum(w_raw)   # probability weights for Σ = E_w[φφ']

# ------------------------------------------------------------
# Pull variables (Missing allowed)
# ------------------------------------------------------------
duration       = hasproperty(df, :duration) ? Vector{Union{Missing,Float64}}(df.duration) : fill(missing, N)
d_wealth_ratio = Vector{Union{Missing,Float64}}(df.d_wealth_ratio)
adj_ratio      = Vector{Union{Missing,Float64}}(df.adj_ratio)
usd_share      = Vector{Union{Missing,Float64}}(df.usd_share)

# -------------------------------
# Common mask across all moments
# -------------------------------
m1 = finmask(duration)
m2 = finmask(d_wealth_ratio)
m4 = finmask(adj_ratio)
m5 = finmask(usd_share)
mask = m1 .& m2 .& m4 .& m5

w2 = w_raw[mask]
w2 ./= sum(w2)
Wobs2 = Diagonal(w2)

dur2 = Float64.(duration[mask])
dwr2 = Float64.(d_wealth_ratio[mask])
adj2 = Float64.(adj_ratio[mask])
usd2 = Float64.(usd_share[mask])

# IFs on common mask (no zero fill)
φ1 = dur2 .- mean(dur2, Weights(w2))

μd = mean(dwr2, Weights(w2))
σd = sum(w2 .* (dwr2 .- μd).^2) / sum(w2)
φ2 = dwr2 .- μd
φ3 = (dwr2 .- μd).^2 .- σd

φ4 = adj2 .- mean(adj2, Weights(w2))

μu = mean(usd2, Weights(w2))
σu = sum(w2 .* (usd2 .- μu).^2) / sum(w2)
φ5 = usd2 .- μu
φ6 = (usd2 .- μu).^2 .- σu

IF2 = hcat(φ1, φ2, φ3, φ4, φ5, φ6)

Σraw = IF2' * Wobs2 * IF2
Σ = Symmetric((Σraw + Σraw')/2)

# -------------------------------
# Prioritized FULL weight matrix
# W = P * (Σ + ridge I)^(-1) * P
# -------------------------------
prio = 2.0
p = ones(6)
p[[1,2,5]] .= prio
P = Diagonal(p)

ridge = 1e-6
Wprio = Symmetric(P * inv(Matrix(Σ) + ridge*I) * P)


# ------------------------------------------------------------
# Save side-by-side (parallel-friendly names)
# ------------------------------------------------------------
isdir(kst.DATA_DIR) || mkpath(kst.DATA_DIR)

writedlm(joinpath(kst.DATA_DIR, "moments_vector_data_weighted_pick6.csv"), datamom)
writedlm(joinpath(kst.DATA_DIR, "moments_names_pick6.csv"), mom_names)

# Σ (covariance) for SEs / reference
writedlm(joinpath(kst.DATA_DIR, "Sigma_pick6.csv"), Matrix(Σ))

# Wprio (weight matrix) for estimation objective
writedlm(joinpath(kst.DATA_DIR, "W_pick6_prio125_diag.csv"), Matrix(Wprio))

println("Saved:")
println("  moments_vector_data_weighted_pick6.csv")
println("  moments_names_pick6.csv")
println("  Sigma_pick6.csv")
println("  W_pick6_prio125_diag.csv  (prio=$(prio) on moments 1,2,5)")
