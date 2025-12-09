using Statistics, LinearAlgebra, DataFrames, CSV
using StatsBase: Weights, mean, var, cov, quantile, std
using Distributions: Normal
using DelimitedFiles: writedlm

include("durable_mod.jl")
include("inflnc_functions.jl")
using Main.kst  

# ---------- read per-observation data ----------
df = CSV.read(joinpath(kst.DATA_DIR, "EFHU_moments_data_weighted.csv"), DataFrame)

N     = nrow(df)
peso  = Vector{Float64}(df.pesoEFHU)
p     = peso ./ sum(peso)  # normalized weights

# Variables (Union{Missing,Float64})
d_value        = Vector{Union{Missing,Float64}}(df.durables)
a_eff          = Vector{Union{Missing,Float64}}(df.a_eff)
usd_share      = Vector{Union{Missing,Float64}}(df.usd_share)
usd_particip   = (.!ismissing.(usd_share)) .& (coalesce.(usd_share, 0.0) .> 0.0)

d_income_ratio = Vector{Union{Missing,Float64}}(df.d_income_ratio)
d_wealth_ratio = Vector{Union{Missing,Float64}}(df.d_wealth_ratio)
adj_ratio      = Vector{Union{Missing,Float64}}(df.adj_ratio)

duration = hasproperty(df, :duration) ?
    Vector{Union{Missing,Float64}}(df.duration) :
    fill(missing, N)

usd_heavy   = Vector{Union{Missing,Float64}}(usd_share .> 0.5)
dwealth_usd = Vector{Union{Missing,Float64}}(ifelse.(usd_particip, d_wealth_ratio, missing))

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

# ---------- influence functions ----------
# Helpers from inflnc_functions.jl:
# aligned_xw, aligned_xyw, wcov, etc.

function IF_wmean(xm, p)
    x, w, mask = aligned_xw(xm, p)
    μ = mean(x, Weights(w))
    out = zeros(length(xm))
    out[mask] = w .* (x .- μ)
    return out
end

function IF_wvar(xm, p)
    x, w, mask = aligned_xw(xm, p)
    μ  = mean(x, Weights(w))
    σ2 = var(x, Weights(w))
    out = zeros(length(xm))
    out[mask] = w .* ((x .- μ).^2 .- σ2)
    return out
end

function IF_wshare_pos(zm, p)
    z, w, mask = aligned_xw(zm, p)
    P = sum(w .* (z .> 0.0))
    out = zeros(length(zm))
    out[mask] = w .* ((z .> 0.0) .- P)
    return out
end

# ---------- Build IFs for the 5 moments ----------
# 1) duration_mean
IF_m1 = IF_wmean(duration, p)

# 2) dwealth_mean
IF_m3 = IF_wmean(d_wealth_ratio, p)

# 3) dwealth_var
IF_m4 = IF_wvar(d_wealth_ratio, p)

# 4) adj_rate (mean adjustment ratio)
IF_m6 = IF_wmean(adj_ratio, p)

# 5) owner_share (d_value > 0)
#IF_m7 = IF_wshare_pos(d_value, p)

#5) Dollar asset shares 
IF_m7 = IF_wmean(usd_share, p)

#6) Dollar Asset Vol 
 IF_m8 = IF_wvar(usd_share, p)


# Stack only these 5 columns
IF_matrix = hcat(
    IF_m1,  # duration_mean
    IF_m3,  # dwealth_mean
    IF_m4,  # dwealth_var
    IF_m6,  # adj_rate
    IF_m7,   # dollar_share
    IF_m8    # dollar_vol

)

# Covariance-style matrix Σ ∝ E[IF IF']
Σ_raw = IF_matrix' * IF_matrix
Σ     = Symmetric((Σ_raw + Σ_raw')/2)

# ---------- save ----------
isdir(kst.DATA_DIR) || mkpath(kst.DATA_DIR)
writedlm(kst.MOMS_FILE, datamom)
writedlm(kst.W_FILE,   Σ)
writedlm(kst.MNAME_FILE, mom_names)

println("✅ Saved 6 moments and Σ ($(size(Σ))). Columns:")
println(join(mom_names, ", "))
