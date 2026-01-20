# ==========================================================================
# 4D MODEL: Moment computation from simulated data
# ==========================================================================

using Statistics

_fin(v) = v[isfinite.(v)]

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
    pd = pea[10]       # durable price
    delta = pea[2]     # depreciation rate

    # Post-burn-in window
    r0 = (sz.burnin - 2):sz.nYears
    r1 = (sz.burnin - 3):(sz.nYears - 1)

    # Extract simulated series
    # In 4D model, we have w (total wealth), s (dollar share), and derived aa, a
    w = simdata.w[r0, :]
    s = simdata.s[r0, :]
    d = simdata.d[r0, :]
    d_lag = simdata.d[r1, :]
    ex = simdata.ex[r0, :]
    y = simdata.y[r0, :]
    adj_r0 = simdata.adjust_indicator[r0, :]
    
    # Derived asset holdings from simulation (or compute if not stored)
    if hasproperty(simdata, :aa) && hasproperty(simdata, :a)
        aa = simdata.aa[r0, :]
        a = simdata.a[r0, :]
    else
        # Derive from w and s
        aa = (1.0 .- s) .* w
        a = (s .* w) ./ max.(ex, 1e-10)
    end
    
    c = hasproperty(simdata, :c) ? simdata.c[r0, :] : zeros(size(w))

    # =========================================================================
    # Moment 1: Duration (years since last adjustment)
    # =========================================================================
    adj_full = simdata.adjust_indicator
    spells_full = running_spells(adj_full)
    
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
    a_loc = aa                       # peso assets
    a_fx = ex .* a                   # dollar assets in peso terms
    a_eff = a_loc .+ a_fx            # total liquid wealth
    
    d_value_vec = vec(Vd)
    a_eff_vec = vec(a_eff)
    
    # Durable-to-total-wealth ratio
    total_wealth = d_value_vec .+ a_eff_vec
    d_wealth_ratio = d_value_vec ./ max.(total_wealth, 1e-12)
    
    dwr_fin = _fin(d_wealth_ratio)
    m2_dwealth_mean = mean(dwr_fin)
    m3_dwealth_var = var(dwr_fin)

    # =========================================================================
    # Moments 5-6: Dollar share of liquid assets
    # =========================================================================
    # Dollar share = dollar assets / total liquid assets
    # In simulation, we track s directly, but can also compute from aa and a
    
    a_fx_vec = vec(a_fx)
    usd_share = a_fx_vec ./ max.(a_eff_vec, 1e-12)
    
    usd_sh_fin = _fin(usd_share)
    m5_dollar_share = mean(usd_sh_fin)
    m6_dollar_vol = var(usd_sh_fin)  # variance, not std

    # =========================================================================
    # Additional moments (for diagnostics, not necessarily matched)
    # =========================================================================
    eps = 1e-6
    c_vec = vec(c)
    cons_vol = std(log.(_fin(c_vec .+ eps)))
    
    d_spend = d .- d_lag .* (1.0 .- delta)
    d_spend_vol = std(_fin(vec(d_spend .+ eps)))
    
    a_eff_vol = std(log.(_fin(a_eff_vec .+ eps)))

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

    # Warn about bad values
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
    end

    return outmoms
end
