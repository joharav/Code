using Statistics, Distributions, LinearAlgebra, DataFrames, CSV, Printf, DelimitedFiles
include("inflnc_functions.jl")


# ---------------------------
# 1. Import the EFHU data
# ---------------------------
dta = CSV.read("../../Stata_Folder/EFHU_moments_data.csv", DataFrame)

# Convert relevant columns to vectors
income         = Vector{Union{Missing, Float64}}(dta.household_income)
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

datamoments = [
    weighted_mean(d_income_ratio, pesoEFHU),
    weighted_mean(d_wealth_ratio, pesoEFHU),
    weighted_mean(adj_ratio, pesoEFHU)
]

println("✅ Moments computed: ", datamoments)
# Filter missing for each variable
valid_income  = .!ismissing.(d_income_ratio)
valid_wealth  = .!ismissing.(d_wealth_ratio)
valid_adj     = .!ismissing.(adj_ratio)

IF_income = zeros(length(d_income_ratio))
IF_income[valid_income] = pesoEFHU[valid_income] .* mean_inflnc(d_income_ratio[valid_income])

IF_wealth = zeros(length(d_wealth_ratio))
IF_wealth[valid_wealth] = pesoEFHU[valid_wealth] .* mean_inflnc(d_wealth_ratio[valid_wealth])

IF_adj = zeros(length(adj_ratio))
IF_adj[valid_adj] = pesoEFHU[valid_adj] .* mean_inflnc(adj_ratio[valid_adj])

IF_matrix = hcat(IF_income, IF_wealth, IF_adj)
@assert size(IF_matrix,1) > 0 && size(IF_matrix,2) > 0 "❌ IF_matrix is empty"

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

writedlm(MOMS_FILE, datamoments)
writedlm(W_FILE, W)

# (Optional) save moment names
mom_names = ["d_income_ratio", "d_wealth_ratio", "adj_ratio"]
writedlm(MNAME_FILE, mom_names)



writedlm("datamoments.txt", datamoments)
writedlm("W_unclustered.txt", W)

# ---------------------------
# 6. Save observation-level data for IFs
# ---------------------------
obs_data = DataFrame(
    hogar_id = hogar_id,
    survey_wave = survey_wave,
    income_q = income_q,
    d_income_ratio = d_income_ratio,
    d_wealth_ratio = d_wealth_ratio,
    adj_ratio = adj_ratio,
    pesoEFHU = pesoEFHU
)

CSV.write(joinpath(DATA_DIR, "EFHU_moments_data_forIF.csv"), obs_data)
println("✅ Exported observation-level data to EFHU_moments_data_forIF.csv")
