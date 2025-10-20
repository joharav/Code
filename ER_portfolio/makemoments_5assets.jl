############  Stata-equivalent moments (m1..m12)  ############
# Matches your Stata code:
# m1  mean(d_income_ratio)
# m2  var(d_income_ratio)
# m3  mean(adj_spell)  where adj_spell = 1/(1+duration_years)
# m4  corr(adj_spell, d_income_ratio)
# m5  corr(adj_spell, usd_sh)
# m6  mean(usd_pt)  with usd_pt = 1(usd_sh>0)
# m7  corr(usd_sh, d_income_ratio)
# m8  corr(usd_sh, a_eff)
# m9  mean(d_wealth_ratio)
# m10 var(d_wealth_ratio)
# m11 mean(duration_years)
# m12 var(log(1 + d_value))

using Statistics

# --- small helpers ---
_fin(v) = v[isfinite.(v)]
function _fin2(x, y)
    m = isfinite.(x) .& isfinite.(y)
    return x[m], y[m]
end

# Time since last adjustment (in periods) for each column at the FINAL row
# adj: T×N indicator (Bool/Int)
function last_spell_lengths(adj::AbstractMatrix)
    T, N = size(adj)
    out = fill(NaN, N)
    @inbounds for j in 1:N
        # walk back from T until the last 1
        len = 0
        for t in T:-1:1
            if adj[t, j] > 0
                break
            end
            len += 1
        end
        out[j] = len
    end
    return out # in periods
end

"""
makemoments_stata(simdata, pea; per_year=4)

Return the 12-moment vector that replicates your Stata definitions.
"""
function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool=false)
    # --- parameters used for valuation ---
    w  = pea[8]   # scale for earnings if needed
    pd = pea[10]  # USD durable price

    # --- aligned windows ---
    t0 = (sz.burnin-2):sz.nYears        # current (T rows)
    t1 = (sz.burnin-3):(sz.nYears-1)    # lagged  (T rows)

    a        = simdata.a[t0, :]
    aa       = simdata.aa[t0, :]
    d        = simdata.d[t0, :]
    ex       = simdata.ex[t0, :]
    c        = simdata.c[t0, :]
    y        = simdata.y[t0, :]
    d_a    = simdata.d_adjust[t0, :]
    adj  = simdata.adjust_indicator[t0, :]
  
    a_lag    = simdata.a[t1, :]
    aa_lag   = simdata.aa[t1, :]
    d_lag    = simdata.d[t1, :]
    ex_lag   = simdata.ex[t1, :]


    # value of durables (local units): price is USD*ex
    Vd   = pd .* ex .* d

    # total liquid savings in local units
    a_loc = aa                      # pesos
    a_fx  = ex .* a                 # dollars valued in local units
    a_eff = a_loc .+ a_fx           # total savings (local units)

    # income in local units (guard against zeros)
    income = max.(w .* y, 1e-12)

    # ---------- Cross-sectional objects (use last period per household) ----------
    # take the final row as the "survey" cross-section
    Vd_cs     = vec(view(Vd, size(Vd,1), :))
    a_eff_cs  = vec(view(a_eff, size(a_eff,1), :))
    income_cs = vec(view(income, size(income,1), :))
    a_fx_cs   = vec(view(a_fx, size(a_fx,1), :))
    a_tot_cs  = max.(a_eff_cs, 1e-12)  # same denom as Stata "total_savings_final"

    # Stata counterparts
    d_value_cs      = Vd_cs
    d_income_ratio  = d_value_cs ./ income_cs
    d_wealth_ratio  = d_value_cs ./ max.(d_value_cs .+ a_eff_cs, 1e-12)
    usd_sh          = a_fx_cs ./ a_tot_cs
    usd_pt          = usd_sh .> 0.0

    # duration (years since last update) at the final date
    last_spells_periods = last_spell_lengths(adj)             # in periods
    duration_years      = last_spells_periods ./ 4
    adj_spell           = 1.0 ./ (1.0 .+ duration_years)      # (0,1]

    # ---------- Moments (exactly as in your Stata block) ----------
    # m1/m2: mean/var of d_income_ratio
    dinc_f = _fin(d_income_ratio)
    m1 = mean(dinc_f)
    m2 = var(dinc_f)

    # m3: mean(adj_spell)
    m3 = mean(_fin(adj_spell))

    # m4: corr(adj_spell, d_income_ratio)
    x, yv = _fin2(adj_spell, d_income_ratio)
    m4 = length(x) > 1 ? cor(x, yv) : NaN

    # m5: corr(adj_spell, usd_sh)
    x, yv = _fin2(adj_spell, usd_sh)
    m5 = length(x) > 1 ? cor(x, yv) : NaN

    # m6: mean(usd_pt)
    m6 = mean(_fin(Float64.(usd_pt)))

    # m7: corr(usd_sh, d_income_ratio)
    x, yv = _fin2(usd_sh, d_income_ratio)
    m7 = length(x) > 1 ? cor(x, yv) : NaN

    # m8: corr(usd_sh, a_eff)
    x, yv = _fin2(usd_sh, a_eff_cs)
    m8 = length(x) > 1 ? cor(x, yv) : NaN

    # m9/m10: mean/var of d_wealth_ratio
    dw_f = _fin(d_wealth_ratio)
    m9  = mean(dw_f)
    m10 = var(dw_f)

    # m11: mean(duration_years)
    m11 = mean(_fin(duration_years))

    # m12: var(log(1 + d_value))
    m12 = var(_fin(log.(1 .+ d_value_cs)))

    # ----- gap–hazard diagnostics -----
    gap_vec, f_x, x_values, h_x, I_d_abs, mu_gap, var_gap, adj_rate_gap =
    adjustment_gaps_sim(d_lag, d_a, adj)

    outmoms = [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12]
    return outmoms, x_values, f_x, h_x
end

