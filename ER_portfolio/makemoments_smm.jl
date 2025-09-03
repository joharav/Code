function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool = false)
  # keep the same front-matter
  outmoms = zeros(sz.nmom)

  beta  = pea[1];  w = pea[8];  pd = pea[10];  rr = 1/beta - 1

  # pull series (both aa and a exist)
  aa               = simdata.aa[sz.burnin-2:sz.nYears, :]
  aa_state         = simdata.aa[sz.burnin-3:sz.nYears-1, :]
  a                = simdata.a[sz.burnin-2:sz.nYears, :]
  a_state          = simdata.a[sz.burnin-3:sz.nYears-1, :]
  d                = simdata.d[sz.burnin-2:sz.nYears, :]
  d_state          = simdata.d[sz.burnin-3:sz.nYears-1, :]
  ex               = simdata.ex[sz.burnin-2:sz.nYears, :]
  c                = simdata.c[sz.burnin-2:sz.nYears, :]
  d_adjust         = simdata.d_adjust[sz.burnin-2:sz.nYears, :]
  adjust_indicator = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]
  y                = simdata.y[sz.burnin-2:sz.nYears, :]

  # effective assets in local currency (for aggregates/ratios)
  a_eff_state  = aa_state .+ ex[1:end-1, :] .* a_state
  a_eff        = aa        .+ ex            .* a

  # adjustment gaps (unchanged)
  adjustment_indicator_vec = vec(adjust_indicator)
  gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio =
      adjustment_gaps_sim(d_state, d_adjust, adjustment_indicator_vec)

  # ===== changes =====
  epsv = 1.0e-12 # guard against 0-denominators

  d_invest = 100 .* (d .- d_state) ./ (d_state .+ epsv)

  # split assets: local aa, foreign (in local units) aF = ex .* a
  aF        = ex .* a
  aF_state  = ex[1:end-1, :] .* a_state

  aa_change = 100 .* (aa .- aa_state) ./ (aa_state .+ epsv)
  aF_change = 100 .* (aF .- aF_state) ./ (aF_state .+ epsv)

  # keep your “effective” asset change too (backwards compatible)
  a_eff_change = 100 .* (a_eff .- a_eff_state) ./ (a_eff_state .+ epsv)

  # consumption change (same)
  c_change = 100 .* (c[2:end, :] .- c[1:end-1, :]) ./ (c[1:end-1, :] .+ epsv)

  # ===== moments =====
  mu_d  = mean(vec(d_invest));  var_d  = var(vec(d_invest))
  mu_c  = mean(vec(c_change));   var_c  = var(vec(c_change))

  # effective asset change (old single-asset proxy)
  mu_a_eff = mean(vec(a_eff_change)); var_a_eff = var(vec(a_eff_change))

  # NEW: local vs foreign (local units) changes
  mu_aa = mean(vec(aa_change));  var_aa = var(vec(aa_change))
  mu_aF = mean(vec(aF_change));  var_aF = var(vec(aF_change))

  mu_d1 = mean(vec(d_state));  var_d1 = var(vec(d_state))
  d_dispersion = var(log1p.(vec(d_state)))

  # durable-value ratios (unchanged)
  ratio_d_income      = vec(pd .* ex .* d) ./ vec(w .* y .+ a_eff_state .* (1 + rr) .+ epsv)
  ratio_d_wealth      = vec(pd .* ex .* d) ./ vec(a_eff_state .* (1 + rr) .+ pd .* ex .* d_state .+ epsv)
  ratio_d_consumption = vec(pd .* ex .* d) ./ vec(c .+ epsv)

  mu_d_income = mean(ratio_d_income)
  mu_d_wealth = mean(ratio_d_wealth)
  mu_d_c      = mean(ratio_d_consumption)

  disp_d_income = compute_dispersion(ratio_d_income)
  disp_d_wealth = compute_dispersion(ratio_d_wealth)
  disp_d_c      = compute_dispersion(ratio_d_consumption)
  disp_d        = compute_dispersion(vec(d))

  IQR_d_income    = disp_d_income[2];  p90_10_d_income = disp_d_income[3]
  IQR_d_wealth    = disp_d_wealth[2];  p90_10_d_wealth = disp_d_wealth[3]
  IQR_d_c         = disp_d_c[2];       p90_10_d_c      = disp_d_c[3]
  IQR_d           = disp_d[2];         p90_10_d        = disp_d[3]

  # optional diagnostics
  corr_d_c  = cor(vec(d), vec(c))
  corr_d_aF = cor(vec(d), vec(aF))
  corr_d_aa = cor(vec(d), vec(aa))

  # A simple local-share moment (how much of financial wealth is local)
  local_share_state = vec( aa_state ./ (a_eff_state .+ epsv) )
  μ_local_share  = mean(local_share_state)
  σ2_local_share = var(local_share_state)

  # plots (unchanged)
  plot_aggregates(simdata)
  d_adjust_time_size(simdata)
  plot_simulated_d_and_a_by_state(simdata)

  # ===== pack moments =====
  # Keep your original first three the same (so old pick=[1,2,3] still works),
  # then append the new aa/aF moments. Fill up to sz.nmom, zeros after.
  core = [
      d_dispersion,            # 1
      mu_d_wealth,             # 2
      adjustment_ratio,        # 3
      mu_aa,                   # 4  (NEW)
      var_aa,                  # 5  (NEW)
      mu_aF,                   # 6  (NEW)
      var_aF,                  # 7  (NEW)
      mu_a_eff,                # 8  (still useful)
      var_a_eff,               # 9
      μ_local_share,           # 10 (NEW)
      σ2_local_share           # 11 (NEW)
  ]
  K = min(length(core), sz.nmom)
  outmoms = zeros(sz.nmom)
  outmoms[1:K] .= core[1:K]

  return outmoms, x_values, f_x, h_x, IQR_d, p90_10_d
end
