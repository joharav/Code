using Statistics, LinearAlgebra, DataFrames, CSV, DelimitedFiles
using StatsBase: Weights, mean, var, cov, quantile
include("inflnc_functions.jl")   # keep if you rely on helpers; we override key IFs below
using Main.kst  # paths: MOMS_FILE, W_FILE, MNAME_FILE, DATA_DIR

# ---------- helpers ----------
_finmask(v) = collect(ismissing.(v) .| .!isfinite.(coalesce.(v, NaN)) .|> (!) )  # finite & not missing

# Weighted KDE density at a point (Silverman bandwidth)
function w_density_at(x::Vector{Float64}, w::Vector{Float64}, θ::Float64)
    n = length(x)
    σ = std(x, Weights(w))
    σ = σ > 0 ? σ : 1e-8
    h = 1.06 * σ * n^(-1/5) + 1e-12
    z = (x .- θ) ./ h
    # Gaussian kernel
    @inbounds return sum(w .* pdf.(Normal(), z)) / h
end

# Weighted correlation components
function wcov(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64})
    μx, μy = mean(x, Weights(w)), mean(y, Weights(w))
    sum(w .* (x .- μx) .* (y .- μy))
end

# ---------- read per-observation data ----------
df = CSV.read(joinpath(kst.DATA_DIR, "EFHU_moments_data_weighted.csv"), DataFrame)

income_q        = Vector{Union{Missing,Float64}}(df.income_q)
income          = copy(income_q)     # for convenience
durables        = Vector{Union{Missing,Float64}}(df.durables)
d_income_ratio  = Vector{Union{Missing,Float64}}(df.d_income_ratio)
d_wealth_ratio  = Vector{Union{Missing,Float64}}(df.d_wealth_ratio)
adj_ratio       = Vector{Union{Missing,Float64}}(df.adj_ratio)
usd_share       = Vector{Union{Missing,Float64}}(df.usd_share)
a_eff           = hasproperty(df, :a_eff) ? Vector{Union{Missing,Float64}}(df.a_eff) : fill(missing, nrow(df))
peso            = Vector{Float64}(df.pesoEFHU)

# Normalize survey weights to probabilities
Wtot = sum(peso)
p = peso ./ Wtot

N = nrow(df)

# Utility: get aligned finite vectors & weights
function aligned_xyw(xu, yu, wu)
    m = [!ismissing(xu[i]) && !ismissing(yu[i]) && isfinite(xu[i]) && isfinite(yu[i]) for i in eachindex(xu)]
    x = Float64.(xu[m]); y = Float64.(yu[m]); w = wu[m]
    return x, y, w, m
end
function aligned_xw(xu, wu)
    m = [!ismissing(xu[i]) && isfinite(xu[i]) for i in eachindex(xu)]
    x = Float64.(xu[m]); w = wu[m]; return x, w, m
end

# ---------- 1) Weighted moments (population) ----------
# m1..m5: d/y distribution
x_dy, w_dy, m_dy = aligned_xw(d_income_ratio, p)
m1  = mean(x_dy, Weights(w_dy))
m2  = var(x_dy, Weights(w_dy))
qs  = quantile(x_dy, [0.25,0.5,0.75]; weights=Weights(w_dy))
m3, m4, m5 = qs

# m6: mean(adj)
x_adj, w_adj, m_adj = aligned_xw(adj_ratio, p)
m6  = mean(x_adj, Weights(w_adj))

# m7: corr(adj, d/y)
x1, x2, w12, m12 = aligned_xyw(adj_ratio, d_income_ratio, p)
μ1, μ2 = mean(x1, Weights(w12)), mean(x2, Weights(w12))
σ1 = sqrt(var(x1, Weights(w12)) + 1e-12)
σ2 = sqrt(var(x2, Weights(w12)) + 1e-12)
m7  = wcov(x1, x2, w12) / (σ1*σ2 + 1e-12)

# m8: corr(adj, usd_share)
xu, xa, wu, mu = aligned_xyw(adj_ratio, usd_share, p)
μx, μu = mean(xu, Weights(wu)), mean(xa, Weights(wu))
σx = sqrt(var(xu, Weights(wu)) + 1e-12)
σu = sqrt(var(xa, Weights(wu)) + 1e-12)
m8  = wcov(xu, xa, wu) / (σx*σu + 1e-12)

# m9: share(usd_share>0)
u_pos = (!ismissing).(usd_share) .& (coalesce.(usd_share, 0.0) .> 0)
m9 = sum(p[u_pos])

# m10: weighted D9–D3 gap by income
xi, wi, mi = aligned_xw(income_q, p)
q30, q90 = quantile(xi, [0.30, 0.90]; weights=Weights(wi))
D3 = (!ismissing).(income_q) .& (Float64.(income_q) .<= q30)
D9 = (!ismissing).(income_q) .& (Float64.(income_q) .>= q90)
yu, wy, my = aligned_xw(usd_share, p)
μ3 = sum(p[D3] .* Float64.(usd_share[D3])) / sum(p[D3])
μ9 = sum(p[D9] .* Float64.(usd_share[D9])) / sum(p[D9])
m10 = μ9 - μ3

# m11: corr(usd_share, d/y)
x3, x4, w34, m34 = aligned_xyw(usd_share, d_income_ratio, p)
μ3x, μ4x = mean(x3, Weights(w34)), mean(x4, Weights(w34))
σ3 = sqrt(var(x3, Weights(w34)) + 1e-12)
σ4 = sqrt(var(x4, Weights(w34)) + 1e-12)
m11 = wcov(x3, x4, w34) / (σ3*σ4 + 1e-12)

# m12: corr(usd_share, a_eff)
if any(.!ismissing.(a_eff))
    x5, x6, w56, m56 = aligned_xyw(usd_share, a_eff, p)
    μ5, μ6 = mean(x5, Weights(w56)), mean(x6, Weights(w56))
    σ5 = sqrt(var(x5, Weights(w56)) + 1e-12)
    σ6 = sqrt(var(x6, Weights(w56)) + 1e-12)
    m12 = wcov(x5, x6, w56) / (σ5*σ6 + 1e-12)
else
    m12 = NaN
end

# m13..m17: d/(d+a)
x_dw, w_dw, m_dw = aligned_xw(d_wealth_ratio, p)
m13 = mean(x_dw, Weights(w_dw))
m14 = var(x_dw, Weights(w_dw))
qdw = quantile(x_dw, [0.25,0.5,0.75]; weights=Weights(w_dw))
m15, m16, m17 = qdw

# m21: spell_mean_y proxy (weighted mean of tenure if present)
tenure = hasproperty(df, :years_since_last_durable_purchase) ? df.years_since_last_durable_purchase :
         (hasproperty(df, :years_since_move) ? df.years_since_move : fill(missing, N))
x_tn, w_tn, m_tn = aligned_xw(tenure, p)
m21 = mean(x_tn, Weights(w_tn))

# m24: Var( log1p(durables) ) weighted
x_d, w_d, m_d = aligned_xw(durables, p)
z = log1p.(x_d)
m24 = var(z, Weights(w_d))

datamoments = [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m21,m24]
mom_names = [
    "d_inc_mean","d_inc_var","d_inc_p25","d_inc_p50","d_inc_p75",
    "adj_rate","corr_adj_dinc","corr_adj_usdsh","usd_particip",
    "usdsh_p90m_p30","corr_usdsh_dinc","corr_usdsh_aeff",
    "d_wealth_mean","d_wealth_var","d_wealth_p25","d_wealth_p50","d_wealth_p75",
    "spell_mean_y","d_dispersion"
]

println("✅ 19 weighted moments computed. length=", length(datamoments))

# ---------- 2) Weighted influence functions (population) ----------
# Probability weights p already normalized (sum p = 1)

# IF for weighted mean: IF_i = p_i * (x_i - μ)
function IF_wmean(x::Vector{Union{Missing,Float64}}, p::Vector{Float64})
    xw, pw, mask = aligned_xw(x, p)
    μ = mean(xw, Weights(pw))
    out = zeros(length(x)); out[mask] = pw .* (xw .- μ); out
end

# IF for weighted variance σ² = Σ p (x-μ)²
function IF_wvar(x::Vector{Union{Missing,Float64}}, p::Vector{Float64})
    xw, pw, mask = aligned_xw(x, p)
    μ = mean(xw, Weights(pw)); σ2 = var(xw, Weights(pw))
    out = zeros(length(x)); out[mask] = pw .* ((xw .- μ).^2 .- σ2); out
end

# IF for weighted quantile θτ: IF_i = p_i*(τ - 1{x≤θ}) / f_w(θ)
function IF_wquant(x::Vector{Union{Missing,Float64}}, p::Vector{Float64}, τ::Float64)
    xw, pw, mask = aligned_xw(x, p)
    θ = quantile(xw, τ; weights=Weights(pw))
    fθ = max(w_density_at(xw, pw, θ), 1e-8)
    infl = pw .* (τ .- (xw .<= θ)) ./ fθ
    out = zeros(length(x)); out[mask] = infl; out
end

# IF for weighted correlation (analogue of unweighted formula with p_i factor)
function IF_wcorr(xm::Vector{Union{Missing,Float64}}, ym::Vector{Union{Missing,Float64}}, p::Vector{Float64})
    x, y, w, mask = aligned_xyw(xm, ym, p)
    μx, μy = mean(x, Weights(w)), mean(y, Weights(w))
    cx, cy = x .- μx, y .- μy
    σx2, σy2 = var(x, Weights(w)), var(y, Weights(w))
    σx = sqrt(σx2 + 1e-12); σy = sqrt(σy2 + 1e-12)
    covxy = wcov(x, y, w)
    ρ = covxy / (σx*σy + 1e-12)
    t1 = (cx .* cy .- covxy) ./ (σx*σy + 1e-12)
    t2 = 0.5 * ρ .* ( (cx.^2 .- σx2) ./ (σx2 + 1e-12) .+ (cy.^2 .- σy2) ./ (σy2 + 1e-12) )
    IFloc = w .* (t1 .- t2)                       # <-- p_i factor
    out = zeros(length(xm)); out[mask] = IFloc; out
end

# IF for weighted participation share P = Σ p * 1{u>0}
function IF_wshare_pos(z::Vector{Union{Missing,Float64}}, p::Vector{Float64})
    zf, pw, mask = aligned_xw(z, p)
    P = sum(pw .* (zf .> 0.0))
    out = zeros(length(z)); out[mask] = pw .* ((zf .> 0.0) .- P); out
end

# IF for weighted D9–D3 with fixed cutpoints
function IF_w_d9d3(y::Vector{Union{Missing,Float64}}, inc::Vector{Union{Missing,Float64}}, p::Vector{Float64})
    xi, wi, mi = aligned_xw(inc, p)
    q30, q90 = quantile(xi, [0.30,0.90]; weights=Weights(wi))
    D3 = (!ismissing).(inc) .& (Float64.(inc) .<= q30)
    D9 = (!ismissing).(inc) .& (Float64.(inc) .>= q90)
    yv, pv, mv = aligned_xw(y, p)
    # map masks; simpler: rebuild on full sample
    μ3 = sum(p[D3] .* Float64.(y[D3])) / sum(p[D3])
    μ9 = sum(p[D9] .* Float64.(y[D9])) / sum(p[D9])
    P3 = sum(p[D3]); P9 = sum(p[D9])
    out = zeros(length(y))
    # contributions only from obs in D9 or D3
    for i in eachindex(y)
        if ismissing(y[i]) || ismissing(inc[i]); continue; end
        if D9[i]
            out[i] = p[i] * (Float64(y[i]) - μ9) / max(P9,1e-12)
        elseif D3[i]
            out[i] = - p[i] * (Float64(y[i]) - μ3) / max(P3,1e-12)
        end
    end
    return out
end

# Build IF columns in moment order
IF_m1  = IF_wmean(d_income_ratio, p)            # mean(d/y)
IF_m2  = IF_wvar(d_income_ratio, p)             # var(d/y)
IF_m3  = IF_wquant(d_income_ratio, p, 0.25)
IF_m4  = IF_wquant(d_income_ratio, p, 0.50)
IF_m5  = IF_wquant(d_income_ratio, p, 0.75)
IF_m6  = IF_wmean(adj_ratio, p)
IF_m7  = IF_wcorr(adj_ratio, d_income_ratio, p)
IF_m8  = IF_wcorr(adj_ratio, usd_share, p)
IF_m9  = IF_wshare_pos(usd_share, p)
IF_m10 = IF_w_d9d3(usd_share, income_q, p)
IF_m11 = IF_wcorr(usd_share, d_income_ratio, p)
IF_m12 = any(.!ismissing.(a_eff)) ? IF_wcorr(usd_share, a_eff, p) : zeros(N)
IF_m13 = IF_wmean(d_wealth_ratio, p)
IF_m14 = IF_wvar(d_wealth_ratio, p)
IF_m15 = IF_wquant(d_wealth_ratio, p, 0.25)
IF_m16 = IF_wquant(d_wealth_ratio, p, 0.50)
IF_m17 = IF_wquant(d_wealth_ratio, p, 0.75)
IF_m21 = IF_wmean(hasproperty(df,:years_since_last_durable_purchase) ? df.years_since_last_durable_purchase :
                   (hasproperty(df,:years_since_move) ? df.years_since_move : fill(missing,N)), p)
IF_m24 = let z = log1p.(Float64.(durables[.!ismissing.(durables)]))
    # reuse IF_wvar on durables through transformation
    idx = findall(.!ismissing.(durables))
    out = zeros(N)
    # construct temporary vector with transformed z at idx
    tmp = Vector{Union{Missing,Float64}}(fill(missing, N)); tmp[idx] .= z
    out = IF_wvar(tmp, p); out
end

IF_matrix = hcat(
    IF_m1, IF_m2, IF_m3, IF_m4, IF_m5,
    IF_m6, IF_m7, IF_m8, IF_m9, IF_m10, IF_m11, IF_m12,
    IF_m13, IF_m14, IF_m15, IF_m16, IF_m17, IF_m21, IF_m24
)
enn = size(IF_matrix,1)

# Weighted Σ: with IFs already scaled by p_i, Σ = sum(IF_i IF_i')
Σ = IF_matrix' * IF_matrix./(enn^2);

# ---------- save ----------
isdir(kst.DATA_DIR) || mkpath(kst.DATA_DIR)
writedlm(kst.MOMS_FILE, datamoments)
writedlm(kst.W_FILE,   Σ)                  # Σ, not inverted
writedlm(kst.MNAME_FILE, mom_names)

println("✅ Saved 19 weighted moments and Σ (size $(size(Σ))).")
