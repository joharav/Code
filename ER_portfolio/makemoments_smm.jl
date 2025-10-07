function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool=false)
  outmoms = zeros(sz.nmom)

  β   = pea[1]; w = pea[8]; pd = pea[10]; θ = pea[16]
  r   = 1/β - 1

  # --- aligned windows ---
  t0 = (sz.burnin-2):sz.nYears        # current (T rows)
  t1 = (sz.burnin-3):(sz.nYears-1)    # lagged  (T rows)
  # T sanity:
  # @assert length(t0) == length(t1)

  # pull series (all T×nFirms)
  a        = simdata.a[t0, :]
  aa       = simdata.aa[t0, :]
  d        = simdata.d[t0, :]
  ex       = simdata.ex[t0, :]
  c        = simdata.c[t0, :]
  y        = simdata.y[t0, :]
  d_adj    = simdata.d_adjust[t0, :]
  adj_ind  = simdata.adjust_indicator[t0, :]

  a_lag    = simdata.a[t1, :]
  aa_lag   = simdata.aa[t1, :]
  d_lag    = simdata.d[t1, :]
  ex_lag   = simdata.ex[t1, :]
  # y_lag  = simdata.y[t1, :]   # only if you need it

  # effective assets (local currency)
  a_eff     = aa .+ ex .* a
  a_eff_lag = aa_lag .+ ex_lag .* a_lag

  # --- adjustment gaps (unchanged) ---
  adj_vec = vec(adj_ind)
  gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio =
      adjustment_gaps_sim(d_lag, d_adj, adj_vec)

  # --- changes (all T×nFirms except c_change which is T-1×nFirms) ---
  d_invest = 100 .* (d .- d_lag) ./ d_lag
  a_change = 100 .* (a_eff .- a_eff_lag) ./ a_eff_lag
  c_change = 100 .* (c[2:end, :] .- c[1:end-1, :]) ./ c[1:end-1, :]

  # --- moments ---
  mu_d  = mean(vec(d_invest));  var_d  = var(vec(d_invest))
  mu_a  = mean(vec(a_change));  var_a  = var(vec(a_change))
  mu_c  = mean(vec(c_change));  var_c  = var(vec(c_change))
  mu_d1 = mean(vec(d_lag));     var_d1 = var(vec(d_lag))
  d_dispersion = var(log1p.(vec(d_lag)))

  # durable ratios
  Vd = pd .* ex .* d
  denom_income = w .* y .+ a_eff_lag .* (1 + r)
  denom_wealth = a_eff_lag .* (1 + r) .+ pd .* ex .* d_lag

  ratio_d_income      = vec(Vd) ./ vec(denom_income)
  ratio_d_wealth      = vec(Vd) ./ vec(denom_wealth)
  ratio_d_consumption = vec(Vd) ./ vec(c)

  mu_d_income = mean(ratio_d_income)
  mu_d_wealth = mean(ratio_d_wealth)
  mu_d_c      = mean(ratio_d_consumption)

  disp_d_income = compute_dispersion(ratio_d_income)
  disp_d_wealth = compute_dispersion(ratio_d_wealth)
  disp_d_c      = compute_dispersion(ratio_d_consumption)
  disp_d        = compute_dispersion(vec(d))

  IQR_d_income    = disp_d_income[2]; p90_10_d_income = disp_d_income[3]
  IQR_d_wealth    = disp_d_wealth[2]; p90_10_d_wealth = disp_d_wealth[3]
  IQR_d_c         = disp_d_c[2];      p90_10_d_c      = disp_d_c[3]
  IQR_d           = disp_d[2];        p90_10_d        = disp_d[3]

  # optional diagnostics
  _corr_dc = cor(vec(d), vec(c))
  _corr_da = cor(vec(d), vec(a_eff))

  # plots (unchanged)
  plot_aggregates(simdata)
  d_adjust_time_size(simdata)
  plot_simulated_d_and_a_by_state(simdata)

  # your selected moments
  outmoms = [d_dispersion, mu_d_wealth, adjustment_ratio]

  return outmoms, x_values, f_x, h_x, IQR_d, p90_10_d
end
