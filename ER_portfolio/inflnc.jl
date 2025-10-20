using Statistics, LinearAlgebra, DataFrames, CSV
using StatsBase: Weights, mean, var, cov, quantile, std
using Distributions: Normal
using DelimitedFiles: writedlm
include("durable_mod.jl")

include("inflnc_functions.jl")
using Main.kst  # expects: DATA_DIR, MOMS_FILE, W_FILE, MNAME_FILE


# ---------- read per-observation data ----------
df = CSV.read(joinpath(kst.DATA_DIR, "EFHU_moments_data_weighted.csv"), DataFrame)

N = nrow(df)
peso  = Vector{Float64}(df.pesoEFHU)
p     = peso ./ sum(peso)  # normalize to probabilities

# Variables (Union{Missing,Float64} where applicable)
d_value        = Vector{Union{Missing,Float64}}(df.durables)
a_eff          = Vector{Union{Missing,Float64}}(df.a_eff)
usd_share      = Vector{Union{Missing,Float64}}(df.usd_share)
usd_particip   = (.!ismissing.(usd_share)) .& (coalesce.(usd_share, 0.0) .> 0.0)  # bool mask for ref
d_income_ratio = Vector{Union{Missing,Float64}}(df.d_income_ratio)
d_wealth_ratio = Vector{Union{Missing,Float64}}(df.d_wealth_ratio)
adj_ratio      = Vector{Union{Missing,Float64}}(df.adj_ratio)
duration       = hasproperty(df, :duration) ? Vector{Union{Missing,Float64}}(df.duration) :
                                             fill(missing, N)

# ---------- data moments (12), computed here OR read the Stata export ----------
# Option A (recommended): take the exact Stata moments to ensure perfect match
mrow = CSV.read(joinpath(kst.DATA_DIR, "moments_vector_data_weighted.csv"), DataFrame)
datamoments = vec(Matrix(mrow)[1, :])  # 12 numbers from your Stata export

# Option B (alt): compute locally (kept here for verification)
function _compute_locally()
    m = Float64[]
    # 1–2: mean/var of d_income_ratio
    x, w, _ = aligned_xw(d_income_ratio, p)
    push!(m, mean(x, Weights(w)))
    push!(m, var(x, Weights(w)))

    # 3: mean(adj_ratio)
    x, w, _ = aligned_xw(adj_ratio, p)
    push!(m, mean(x, Weights(w)))

    # 4: corr(adj_ratio, d_income_ratio)
    x1, x2, w12, _ = aligned_xyw(adj_ratio, d_income_ratio, p)
    σ1 = sqrt(var(x1, Weights(w12)) + 1e-12)
    σ2 = sqrt(var(x2, Weights(w12)) + 1e-12)
    ρ  = wcov(x1, x2, w12) / (σ1 * σ2 + 1e-12)
    push!(m, ρ)

    # 5: corr(adj_ratio, usd_share)
    x1, x2, w12, _ = aligned_xyw(adj_ratio, usd_share, p)
    σ1 = sqrt(var(x1, Weights(w12)) + 1e-12)
    σ2 = sqrt(var(x2, Weights(w12)) + 1e-12)
    push!(m, wcov(x1, x2, w12) / (σ1 * σ2 + 1e-12))

    # 6: share(usd_share>0)
    mask = (!ismissing).(usd_share) .& (coalesce.(usd_share, 0.0) .> 0.0)
    push!(m, sum(p[mask]))

    # 7: corr(usd_share, d_income_ratio)
    x1, x2, w12, _ = aligned_xyw(usd_share, d_income_ratio, p)
    σ1 = sqrt(var(x1, Weights(w12)) + 1e-12)
    σ2 = sqrt(var(x2, Weights(w12)) + 1e-12)
    push!(m, wcov(x1, x2, w12) / (σ1 * σ2 + 1e-12))

    # 8: corr(usd_share, a_eff)
    if any(.!ismissing.(a_eff))
        x1, x2, w12, _ = aligned_xyw(usd_share, a_eff, p)
        σ1 = sqrt(var(x1, Weights(w12)) + 1e-12)
        σ2 = sqrt(var(x2, Weights(w12)) + 1e-12)
        push!(m, wcov(x1, x2, w12) / (σ1 * σ2 + 1e-12))
    else
        push!(m, NaN)
    end

    # 9–10: mean/var of d_wealth_ratio
    x, w, _ = aligned_xw(d_wealth_ratio, p)
    push!(m, mean(x, Weights(w)))
    push!(m, var(x, Weights(w)))

    # 11: mean(duration)
    x, w, _ = aligned_xw(duration, p)
    push!(m, mean(x, Weights(w)))

    # 12: var(log1p(d_value))
    # (transform then weighted variance)
    x, w, _ = aligned_xw(d_value, p)
    z = log1p.(x)
    push!(m, var(z, Weights(w)))

    return m
end
# localcheck = _compute_locally()  # uncomment to verify vs Stata’s `datamoments`

mom_names = [
    "d_inc_mean","d_inc_var",
    "adj_rate",
    "corr_adj_dinc","corr_adj_usdsh",
    "usd_particip",
    "corr_usdsh_dinc","corr_usdsh_aeff",
    "d_wealth_mean","d_wealth_var",
    "spell_mean_y",
    "d_log1p_var"
]

# ---------- influence functions for the 12 moments ----------
# IF for weighted mean: IF_i = p_i*(x_i - μ)
function IF_wmean(xm, p::Vector{Float64})
    x, w, mask = aligned_xw(xm, p)
    μ = mean(x, Weights(w))
    out = zeros(length(xm)); out[mask] = w .* (x .- μ); out
end

# IF for weighted variance (population): σ² = Σ p (x-μ)²
function IF_wvar(xm, p::Vector{Float64})
    x, w, mask = aligned_xw(xm, p)
    μ  = mean(x, Weights(w))
    σ2 = var(x, Weights(w))
    out = zeros(length(xm)); out[mask] = w .* ((x .- μ).^2 .- σ2); out
end

# IF for weighted correlation ρ_xy
function IF_wcorr(xm, ym, p::Vector{Float64})
    x, y, w, mask = aligned_xyw(xm, ym, p)
    μx, μy = mean(x, Weights(w)), mean(y, Weights(w))
    cx, cy = x .- μx, y .- μy
    σx2, σy2 = var(x, Weights(w)), var(y, Weights(w))
    σx = sqrt(σx2 + 1e-12); σy = sqrt(σy2 + 1e-12)
    covxy = wcov(x, y, w)
    ρ = covxy / (σx*σy + 1e-12)
    t1 = (cx .* cy .- covxy) ./ (σx*σy + 1e-12)
    t2 = 0.5 * ρ .* ((cx.^2 .- σx2) ./ (σx2 + 1e-12) .+ (cy.^2 .- σy2) ./ (σy2 + 1e-12))
    loc = w .* (t1 .- t2)
    out = zeros(length(xm)); out[mask] = loc; out
end

# IF for weighted participation share P = Σ p * 1{u>0}
function IF_wshare_pos(zm, p::Vector{Float64})
    z, w, mask = aligned_xw(zm, p)
    P = sum(w .* (z .> 0.0))
    out = zeros(length(zm)); out[mask] = w .* ((z .> 0.0) .- P); out
end

# Build IF columns in the exact order of mom_names (12 cols)
IF_m1  = IF_wmean(d_income_ratio, p)                     # d_inc_mean
IF_m2  = IF_wvar(d_income_ratio, p)                      # d_inc_var
IF_m3  = IF_wmean(adj_ratio, p)                          # adj_rate
IF_m4  = IF_wcorr(adj_ratio, d_income_ratio, p)          # corr_adj_dinc
IF_m5  = IF_wcorr(adj_ratio, usd_share, p)               # corr_adj_usdsh
IF_m6  = IF_wshare_pos(usd_share, p)                     # usd_particip
IF_m7  = IF_wcorr(usd_share, d_income_ratio, p)          # corr_usdsh_dinc
IF_m8  = any(.!ismissing.(a_eff)) ? IF_wcorr(usd_share, a_eff, p) : zeros(N)  # corr_usdsh_aeff
IF_m9  = IF_wmean(d_wealth_ratio, p)                     # d_wealth_mean
IF_m10 = IF_wvar(d_wealth_ratio, p)                      # d_wealth_var
IF_m11 = IF_wmean(duration, p)                           # spell_mean_y
# m12: var(log1p(d_value))
x, _, idx = aligned_xw(d_value, p)
z = log1p.(x)
tmp = Vector{Union{Missing,Float64}}(fill(missing, N))
tmp[idx] = z
IF_m12 = IF_wvar(tmp, p)


IF_matrix = hcat(
    IF_m1, IF_m2, IF_m3, IF_m4, IF_m5, IF_m6,
    IF_m7, IF_m8, IF_m9, IF_m10, IF_m11, IF_m12
)

# ---------- covariance of moments ----------
# IMPORTANT scaling note:
# - IF columns above already include the probability weight p_i.
# - The (population) variance of the moment vector m = Σ_i IF_i is Σ_i IF_i IF_i'.
# - If you prefer "sample-average" scaling, divide by N^2 (your prior code).
#Σ = IF_matrix' * IF_matrix               # population-style
 Σ = (IF_matrix' * IF_matrix) / (N^2)   # sample-average style (matches your earlier line)

# ---------- save ----------
isdir(kst.DATA_DIR) || mkpath(kst.DATA_DIR)
writedlm(kst.MOMS_FILE, datamoments)
writedlm(kst.W_FILE,   Σ)
writedlm(kst.MNAME_FILE, mom_names)

println("✅ Saved 12 moments and Σ ($(size(Σ))). Columns:")
println(join(mom_names, ", "))
