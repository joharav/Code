function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool = false)

# Initialize the output moments vector
outmoms = zeros(sz.nmom)

# Constants from `pea`
beta  = pea[1]
w     = pea[8]
pd    = pea[10]
theta = pea[16]
rr = 1/beta - 1

# Extract variables from simulation data
a                = simdata.a[sz.burnin-2:sz.nYears, :]
a_state          = simdata.a[sz.burnin-3:sz.nYears-1, :]
d                = simdata.d[sz.burnin-2:sz.nYears, :]
d_state          = simdata.d[sz.burnin-3:sz.nYears-1, :]
ex               = simdata.ex[sz.burnin-2:sz.nYears, :]
c                = simdata.c[sz.burnin-2:sz.nYears, :]
d_adjust         = simdata.d_adjust[sz.burnin-2:sz.nYears, :]
adjust_indicator = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]
y                = simdata.y[sz.burnin-2:sz.nYears, :]
a_effective      = theta .* ex .* a_state .+ (1 - theta) .* a_state

# Calculate adjustment gaps
adjustment_indicator = vec(adjust_indicator)
gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio =
adjustment_gaps_sim(d_state, d_adjust, adjustment_indicator)

# Changes
d_invest = 100*(d .- d_state) ./ d_state
a_change = 100*(a .- a_state) ./ a_state
c_change = 100*(c[2:end] .- c[1:end-1]) ./ c[1:end-1]

# Moments calculations
mu_d  = mean(vec(d_invest))
var_d = var(vec(d_invest))
mu_a  = mean(vec(a_change))
var_a = var(vec(a_change))
mu_c  = mean(vec(c_change))
var_c = var(vec(c_change))
mu_d1 = mean(vec(d_state))
var_d1 = var(vec(d_state))

# Ratios for durable share
ratio_d_income      = vec(pd .* ex .* d) ./ vec(w .* y .+ a_effective .* (1 + rr))
ratio_d_wealth      = vec(pd .* ex .* d) ./ vec(a_effective .* (1 + rr) .+ pd .* ex .* d_state)
ratio_d_consumption = vec(pd .* ex .* d) ./ vec(c)

mu_d_income = mean(ratio_d_income)
mu_d_wealth = mean(ratio_d_wealth)
mu_d_c      = mean(ratio_d_consumption)

# Dispersion measures
disp_d_income = compute_dispersion(ratio_d_income)
disp_d_wealth = compute_dispersion(ratio_d_wealth)
disp_d_c      = compute_dispersion(ratio_d_consumption)
disp_d        = compute_dispersion(vec(d))

# Extract IQR and P90/P10
IQR_d_income      = disp_d_income[2]
IQR_d_wealth      = disp_d_wealth[2]
IQR_d_c           = disp_d_c[2]
IQR_d             = disp_d[2]

p90_10_d_income   = disp_d_income[3]
p90_10_d_wealth   = disp_d_wealth[3]
p90_10_d_c        = disp_d_c[3]
p90_10_d          = disp_d[3]

# Optional: correlations for shock identification
corr_d_c = cor(vec(d), vec(c))
corr_d_a = cor(vec(d), vec(a))


# # Populate outmoms
# outmoms[1]  = mu_d
# outmoms[2]  = var_d
# outmoms[3]  = mu_a
# outmoms[4]  = var_a
# outmoms[5]  = mu_c
# outmoms[6]  = var_c
# outmoms[7]  = mu_d_income
# outmoms[8]  = mu_d_wealth
# outmoms[9]  = mu_d_c
# outmoms[10] = mu_gap
# outmoms[11] = var_gap
# outmoms[12] = adjustment_ratio
# outmoms[13] = IQR_d_income
# outmoms[14] = IQR_d_wealth
# outmoms[15] = IQR_d_c
# outmoms[16] = p90_10_d_income
# outmoms[17] = p90_10_d_wealth
# outmoms[18] = p90_10_d_c
# outmoms[19] = corr_d_c
# outmoms[20] = corr_d_a
# outmoms[21] = cev  # optional if you compute welfare

outmoms = [mu_d_income, mu_d_wealth, adjustment_ratio]  # Select moments based on sz.pick


return outmoms, x_values, f_x, h_x, IQR_d, p90_10_d
end
