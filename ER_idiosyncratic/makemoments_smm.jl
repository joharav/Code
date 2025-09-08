function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool = false)

# Initialize the output moments vector
outmoms = zeros(sz.nmom)

# Constants from `pea`
beta  = pea[1]
w     = pea[8]
pd    = pea[10]
theta = pea[16]
rr = 1/beta - 1
delta = pea[2]

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

# durable flow x_t = d_t - (1-δ)d_{t-1}
x_flow = d[2:end, :] .- (1 - delta) .* d[1:end-1, :]

# NEW: volatilities to use in the counterfactual table
vol_c = std(vec(c_change))
vol_x = std(vec(x_flow))

# Moments calculations
mu_d  = mean(vec(d_invest))
var_d = var(vec(d_invest))
mu_a  = mean(vec(a_change))
var_a = var(vec(a_change))
mu_c  = mean(vec(c_change))
var_c = var(vec(c_change))
mu_d1 = mean(vec(d_state))
var_d1 = var(vec(d_state))
d_dispersion = var(log1p.(vec(d_state)))

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

plot_aggregates(simdata)
d_adjust_time_size(simdata)
plot_simulated_d_and_a_by_state(simdata)
if settings.verbose==true

  #  plotgaps(x_values, f_x, h_x, gap_vec; shock=shock)
#    plotdensities(x_values_d_income, f_d_income, "f_income"; shock=shock)
#    plotdensities(x_values_d_wealth, f_d_wealth, "f_wealth"; shock=shock)
#    plotdensities(x_values_d_consumption, f_d_consumption, "d_c"; shock=shock)
    plot_aggregates(simdata)
  #  d_adjust_time_size(simdata)
  #  plot_simulated_d_and_a_by_state(simdata)
end

outmoms = [d_dispersion, mu_d_wealth, adjustment_ratio,mu_a, vol_c, vol_x]  # Select moments based on sz.pick


return outmoms, x_values, f_x, h_x, IQR_d, p90_10_d
end
