# ==========================================================================
# 4D MODEL: Moment computation from simulated data
# ==========================================================================

using Statistics
using Printf

_fin(v) = v[isfinite.(v)]
_pos(v) = v[isfinite.(v) .& (v .> 0.0)]   # finite AND strictly positive

# Running time-since-last-adjustment (in periods) for each (t,i)
function running_spells(adj::AbstractMatrix)
    T, N = size(adj)
    out = zeros(Float64, T, N)
    @inbounds for j in 1:N
        s = 0.0
        for t in 1:T
            if adj[t, j] > 0
                s = 0.0
            else
                s += 1.0
            end
            out[t, j] = s
        end
    end
    return out
end


function makemoments(simdata::NamedTuple, pea::Vector{Float64};
                     per_year::Int = 4, shock::Bool = false)

    # Parameters
    pd    = pea[10]      # durable price
    delta = pea[2]       # depreciation rate (used only for diagnostics below)

    # Post-burn-in window
    r0 = (sz.burnin - 2):sz.nYears
    r1 = (sz.burnin - 3):(sz.nYears - 1)

    # Extract simulated series
    w      = simdata.w[r0, :]
    s      = simdata.s[r0, :]
    d      = simdata.d[r0, :]
    d_lag  = simdata.d[r1, :]
    ex     = simdata.ex[r0, :]
    adj_r0 = simdata.adjust_indicator[r0, :]
    d_adjust= simdata.d_adjust[r0, :]

    # Consumption (optional; only used for diagnostics)
    c = hasproperty(simdata, :c) ? simdata.c[r0, :] : zeros(size(w))

    # IMPORTANT: derive aa and a from (w,s,ex) for consistency
    # (do not trust stored aa/a unless you're 100% sure they are clamped & consistent)
    aa = (1.0 .- s) .* w
    a  = (s .* w) ./ max.(ex, 1e-12)

    # =========================================================================
    # Moment 1: Duration (years since last adjustment)
    # =========================================================================
    spells_full = running_spells(simdata.adjust_indicator)
    t_survey = last(r0)
    duration_years_cs = spells_full[t_survey, :] ./ per_year
    m1_duration_mean = mean(_fin(vec(duration_years_cs)))

    # =========================================================================
    # Moment 4: Adjustment rate (annualized)
    # =========================================================================
    q_q = mean(vec(adj_r0 .> 0.0))                # quarterly probability
    m4_adj_rate = 1.0 - (1.0 - q_q)^per_year      # annual probability

    # =========================================================================
    # Moments 2-3: Durable wealth ratio
    # =========================================================================
    # Durable value in local currency
    Vd = pd .* ex .* d

    # Total liquid wealth in local currency
    a_loc = aa                   # peso assets
    a_fx  = ex .* a              # dollar assets in peso terms
    a_eff = a_loc .+ a_fx        # total liquid wealth (peso)

    d_value_vec = vec(Vd)
    a_eff_vec   = vec(a_eff)

    total_wealth   = d_value_vec .+ a_eff_vec
    d_wealth_ratio = d_value_vec ./ max.(total_wealth, 1e-12)

    dwr_fin = _fin(d_wealth_ratio)
    m2_dwealth_mean = mean(dwr_fin)
    m3_dwealth_var  = var(dwr_fin)

    # =========================================================================
    # Moments 5-6: Dollar share of liquid assets
    # =========================================================================
    a_fx_vec   = vec(a_fx)
    usd_share  = a_fx_vec ./ max.(a_eff_vec, 1e-12)

    usd_sh_fin = _fin(usd_share)
    m5_dollar_share = mean(usd_sh_fin)
    m6_dollar_vol   = var(usd_sh_fin)   # variance (not std)

    # =========================================================================
    # Additional diagnostics (never crash on log of non-positive)
    # =========================================================================
    # If you don't need these, delete this whole block.
    c_vec = vec(c)

    c_pos = _pos(c_vec)
    cons_vol = isempty(c_pos) ? NaN : std(log.(c_pos))

    d_spend = d .- d_lag .* (1.0 .- delta)
    d_spend_vol = std(_fin(vec(d_spend)))

    a_pos = _pos(a_eff_vec)
    a_eff_vol = isempty(a_pos) ? NaN : std(log.(a_pos))

    # =========================================================================
    # Output vector (matches data moment order)
    # Order: duration_mean, dwealth_mean, dwealth_var, adj_rate, dollar_share, dollar_vol
    # =========================================================================
    outmoms = [
        m1_duration_mean,   # 1: duration_mean (years)
        m2_dwealth_mean,    # 2: dwealth_mean
        m3_dwealth_var,     # 3: dwealth_var
        m4_adj_rate,        # 4: adj_rate (annualized)
        m5_dollar_share,    # 5: dollar_share
        m6_dollar_vol,      # 6: dollar_vol
    ]

    if any(.!isfinite.(outmoms))
        badix = findfirst(!isfinite, outmoms)
        @warn "makemoments: NaN/Inf detected" moment=badix value=outmoms[badix]
    end

    if settings.verbose
        println("\n=== Simulated Moments ===")
        mom_names = ["duration_mean", "dwealth_mean", "dwealth_var",
                     "adj_rate", "dollar_share", "dollar_vol"]
        for (i, (name, val)) in enumerate(zip(mom_names, outmoms))
            @printf("  %2d. %-15s = %.6f\n", i, name, val)
        end

        println("\n=== Diagnostics (not matched) ===")
        @printf("  cons_vol(log c)      = %s\n", string(cons_vol))
        @printf("  d_spend_vol          = %.6f\n", d_spend_vol)
        @printf("  a_eff_vol(log a_eff) = %s\n", string(a_eff_vol))
        @printf("  min(w)=%g  min(d)=%g  min(ex)=%g  min(s)=%g\n",
                minimum(w), minimum(d), minimum(ex), minimum(s))
        @printf("  min(a_eff)=%g  min(c)=%g\n", minimum(a_eff), minimum(c))
    end


    gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio =adjustment_gaps_sim(d_lag,d_adjust,adj_r0)



    return outmoms, x_values, f_x, h_x, gap_vec
end
