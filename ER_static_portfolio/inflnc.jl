using Statistics
using StatsBase: Weights, mean, var

# Align x and weights, drop missing/nonfinite
function aligned_xw(xu::AbstractVector, w::Vector{Float64})
    mask = [!ismissing(xu[i]) && isfinite(Float64(xu[i])) for i in eachindex(xu)]
    x  = Float64.(xu[mask])
    ww = w[mask]          # DO NOT renormalize here
    return x, ww, mask
end


"Weighted mean IF: φ_i = x_i - μ_w  (NOT multiplied by weights)"
function IF_wmean(xu::AbstractVector, w::Vector{Float64})
    x, ww, mask = aligned_xw(xu, w)
    μ = mean(x, Weights(ww))
    out = zeros(Float64, length(xu))
    out[mask] = x .- μ
    return out
end

"Weighted variance IF (prob-weighted, no finite-sample correction)"
function IF_wvar(xu::AbstractVector, w::Vector{Float64})
    x, ww, mask = aligned_xw(xu, w)
    μ  = mean(x, Weights(ww))
    σ2 = sum(ww .* (x .- μ).^2) / sum(ww)   # FIX
    out = zeros(Float64, length(xu))
    out[mask] = (x .- μ).^2 .- σ2
    return out
end

function IF_wshare_pos(xu::AbstractVector, w::Vector{Float64})
    x, ww, mask = aligned_xw(xu, w)
    ind = x .> 0.0
    θ = sum(ww .* ind) / sum(ww)            # FIX
    out = zeros(Float64, length(xu))
    out[mask] = Float64.(ind) .- θ
    return out
end

using LinearAlgebra
using DataFrames, CSV
using DelimitedFiles: writedlm

include("durable_mod.jl")
using Main.kst

df = CSV.read(joinpath(kst.DATA_DIR, "EFHU_moments_data_weighted.csv"), DataFrame)



# ---------- read data moment vector ----------
mrow_full = CSV.read(joinpath(kst.DATA_DIR, "moments_vector_data_weighted.csv"), DataFrame)
datamom_full = vec(Matrix(mrow_full)[1, :])  # length 10 in original order

# Original 10-moment order (for reference):
# 1: duration_mean
# 2: usd_particip_mean
# 3: dwealth_mean
# 4: dwealth_var
# 5: usd_particip_var
# 6: adj_rate
# 7: owner_share
# 8: usd_share_mean
# 9: usd_heavy_share
# 10: dwealth_cond_usd

# We now only want 5 moments: 1,3,4,6,7
const pick6 = [1, 3, 4, 6, 8, 5]
datamom = datamom_full[pick6]

mom_names = [
    "duration_mean",   # m1
    "dwealth_mean",    # m3
    "dwealth_var",     # m4
    "adj_rate",        # m6
    "dollar_share",      # m7
    "dollar_vol"        # m8
]


N = nrow(df)
w_raw = Float64.(df.pesoEFHU)
w = w_raw ./ sum(w_raw)   # probability weights

# Pull variables (allow Missing)
duration       = hasproperty(df, :duration) ? Vector{Union{Missing,Float64}}(df.duration) : fill(missing, N)
d_wealth_ratio = Vector{Union{Missing,Float64}}(df.d_wealth_ratio)
adj_ratio      = Vector{Union{Missing,Float64}}(df.adj_ratio)
usd_share      = Vector{Union{Missing,Float64}}(df.usd_share)

# IF columns (6 moments)
IF_m1 = IF_wmean(duration, w)           # duration_mean
IF_m2 = IF_wmean(d_wealth_ratio, w)     # dwealth_mean
IF_m3 = IF_wvar(d_wealth_ratio, w)      # dwealth_var
IF_m4 = IF_wmean(adj_ratio, w)          # adj_rate (must match how you define it in data moments)
IF_m5 = IF_wmean(usd_share, w)          # dollar_share
IF_m6 = IF_wvar(usd_share, w)           # dollar_vol

IF = hcat(IF_m1, IF_m2, IF_m3, IF_m4, IF_m5, IF_m6)

# Weighted covariance of moments: Σ = E_w[φ φ']
# Implement as φ' diag(w) φ with w normalized to sum to 1 over ALL obs,
# but note: each column has its own missing mask. That’s fine; zeros mean dropped obs.
W = Diagonal(w)
Σ = Symmetric((IF' * W * IF + (IF' * W * IF)') / 2)





# Save
isdir(kst.DATA_DIR) || mkpath(kst.DATA_DIR)
writedlm(kst.MOMS_FILE, datamom)
writedlm(kst.W_FILE, Matrix(Σ))
writedlm(kst.MNAME_FILE, mom_names)

println("Saved 6 moments and Σ = $(size(Σ)) with names: ", join(mom_names, ", "))
