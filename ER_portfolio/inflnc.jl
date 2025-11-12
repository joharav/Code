using Statistics, LinearAlgebra, DataFrames, CSV
using StatsBase: Weights, mean, var, cov, quantile, std
using Distributions: Normal
using DelimitedFiles: writedlm
include("durable_mod.jl")
include("inflnc_functions.jl")
using Main.kst  
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
datamoments = vec(Matrix(mrow)[1, :])  
const pick6= [11,6, 9, 10,12]  # moments to use (1-based)
datamoments = datamoments[pick6]  # select only the 6 moments we use

mom_names = [
    "spell_mean_y",
    "usd_particip",
    "d_wealth_mean","d_wealth_var",
    "usd_particip_var"]

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

IF_m1  = IF_wmean(duration, p)                          # duration
IF_m2  = IF_wshare_pos(usd_share, p)                     # usd_particip
IF_m3  = IF_wmean(d_wealth_ratio, p)                     # d_wealth_mean
IF_m4 = IF_wvar(d_wealth_ratio, p)                      # d_wealth_var
IF_m5  = IF_wvar(usd_share, p)          # usd_particip_var


IF_matrix = hcat(
    IF_m1, IF_m2, IF_m3, IF_m4, IF_m5
)

#Σ = IF_matrix' * IF_matrix               # population-style
 Σ = (IF_matrix' * IF_matrix)  # sample-average style (matches your earlier line)
 Σ = Symmetric((Σ+Σ')/2)  # ensure symmetry

# ---------- save ----------
isdir(kst.DATA_DIR) || mkpath(kst.DATA_DIR)
writedlm(kst.MOMS_FILE, datamoments)
writedlm(kst.W_FILE,   Σ)
writedlm(kst.MNAME_FILE, mom_names)

println("✅ Saved 5 moments and Σ ($(size(Σ))). Columns:")
println(join(mom_names, ", "))
