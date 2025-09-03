using Statistics, Distributions, LinearAlgebra, DataFrames, CSV, Printf, DelimitedFiles
include("inflnc_functions.jl")
using Main.kst


# ---------------------------
# 1. Import the EFHU data
# ---------------------------
# Current working directory is assumed
dta = CSV.read("Data/EFHU_moments_data.csv", DataFrame)

# Convert relevant columns to vectors
income         = Vector{Union{Missing, Float64}}(dta.household_income_q)
durables       = Vector{Union{Missing, Float64}}(dta.durables_final)
d_income_ratio = Vector{Union{Missing, Float64}}(dta.d_income_ratio)
d_wealth_ratio = Vector{Union{Missing, Float64}}(dta.d_wealth_ratio)
adj_ratio      = Vector{Union{Missing, Float64}}(dta.adj_ratio)
pesoEFHU       = Vector{Float64}(dta.pesoEFHU)  # assuming no missing here
income_q      = Vector{Union{Missing, Float64}}(dta.income_q)
survey_wave       = Vector{Float64}(dta.survey_wave)  # assuming no missing here
hogar_id      = Vector{Union{Missing, Float64}}(dta.hogar_id)

# ---------------------------
# 2. Compute weighted moments
# ---------------------------
function weighted_mean(x::Vector{Union{Missing,Float64}}, w::Vector{Float64})
    valid = .!ismissing.(x)
    return sum(x[valid] .* w[valid]) / sum(w[valid])
end

function weighted_var(x::Vector{Union{Missing,Float64}}, w::Vector{Float64})
    v  = .!ismissing.(x)
    xv = Float64.(x[v])   # <- force Float64
    wv = w[v]
    μ  = sum(wv .* xv) / sum(wv)
    return sum(wv .* (xv .- μ).^2) / sum(wv)
end

# ---------- d_dispersion = log(1 + durables) ----------
d_dispersion = Vector{Union{Missing,Float64}}(undef, length(durables))
fill!(d_dispersion, missing)  # <- important to avoid undefined reads
valid_d = .!ismissing.(durables) .& (coalesce.(durables, 0.0) .>= 0)
d_dispersion[valid_d] = log1p.(Float64.(durables[valid_d]))


datamoments = [
    weighted_var(d_dispersion, pesoEFHU),
    weighted_mean(d_wealth_ratio, pesoEFHU),
    weighted_mean(adj_ratio, pesoEFHU)
]

println("✅ Moments computed: ", datamoments)
## ----- Influence functions (use masks; cast to Float64) -----
valid_dispersion = .!ismissing.(d_dispersion)
valid_wealth     = .!ismissing.(d_wealth_ratio)
valid_adj        = .!ismissing.(adj_ratio)

IF_dispersion = zeros(length(d_dispersion))
x_disp = Float64.(d_dispersion[valid_dispersion])
IF_dispersion[valid_dispersion] =
    pesoEFHU[valid_dispersion] .* var_inflnc(x_disp)

IF_wealth = zeros(length(d_wealth_ratio))
x_wealth = Float64.(d_wealth_ratio[valid_wealth])
IF_wealth[valid_wealth] =
    pesoEFHU[valid_wealth] .* mean_inflnc(x_wealth)

IF_adj = zeros(length(adj_ratio))
x_adj = Float64.(adj_ratio[valid_adj])
IF_adj[valid_adj] =
    pesoEFHU[valid_adj] .* mean_inflnc(x_adj)

IF_matrix = hcat(IF_dispersion, IF_wealth, IF_adj)

println("✅ Influence function matrix size: ", size(IF_matrix))
# ---------------------------
# 4. Weighting matrix for SMM
# ---------------------------
nobs = size(IF_matrix,1)
W = IF_matrix' * IF_matrix / nobs^2       # simple (unclustered) weighting
@assert size(W,1) > 0 && size(W,2) > 0 "❌ Weighting matrix W is empty"

println("✅ Weighting matrix size: ", size(W))

# ---------------------------
# 5. Save matrices for SMM
# ---------------------------

# Make sure folders exist
isdir("Data") || mkdir("Data")
isdir("Output") || mkdir("Output")

writedlm(kst.MOMS_FILE, datamoments)
writedlm(kst.W_FILE, W)

# (Optional) save moment names
mom_names = ["d_dispersion", "d_wealth_ratio", "adj_ratio"]
writedlm(kst.MNAME_FILE, mom_names)



writedlm("datamoments.txt", datamoments)
writedlm("W_unclustered.txt", W)

# ---------------------------
# 6. Save observation-level data for IFs
# ---------------------------
obs_data = DataFrame(
    hogar_id = hogar_id,
    survey_wave = survey_wave,
    income_q = income_q,
    income = income,
    d_dispersion = d_dispersion,
    d_wealth_ratio = d_wealth_ratio,
    adj_ratio = adj_ratio,
    pesoEFHU = pesoEFHU
)

CSV.write(joinpath(kst.DATA_DIR, "EFHU_moments_data_forIF.csv"), obs_data)
println("✅ Exported observation-level data to EFHU_moments_data_forIF.csv")
