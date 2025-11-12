using Statistics

# ---------- helpers ----------
_fin(v) = v[isfinite.(v)]
@inline _safe_den(x) = max(x, 1e-12)

@inline function safe_cor(x::AbstractVector, y::AbstractVector)
    # coerce Bool→Float64 and drop non-finite pairs
    xx = Float64.(x); yy = Float64.(y)
    m = isfinite.(xx) .& isfinite.(yy)
    if !any(m) || std(xx[m]) ≤ 0 || std(yy[m]) ≤ 0
        return 0.0    # or NaN if you prefer to flag
    end
    return cor(xx[m], yy[m])
end


# Running time-since-last-adjustment (in periods) for each (t,i) in r0
# adj: T×N Bool/Int in the selected window r0
function running_spells(adj::AbstractMatrix)
    T, N = size(adj)
    out = zeros(Float64, T, N)
    @inbounds for j in 1:N
        s = 0.0
        for t in 1:T
            if adj[t,j] > 0
                s = 0.0            # reset on adjustment
            else
                s += 1.0
            end
            out[t,j] = s
        end
    end
    return out  # same size as adj
end


function makemoments(simdata::NamedTuple, pea::Vector{Float64};
                     per_year::Int = 4, shock::Bool=false)

    # --- parameters ---
    w  = pea[8]    # earnings scale
    pd = pea[10]   # durable price (USD)

    # --- post-burn-in window ---
    r0 = (sz.burnin-2):sz.nYears
    r1 = (sz.burnin-3):(sz.nYears-1)

    # slice panel
    a      = simdata.a[r0, :]
    aa     = simdata.aa[r0, :]
    d      = simdata.d[r0, :]
    ex     = simdata.ex[r0, :]
    y      = simdata.y[r0, :]
#    c      = simdata.c[r0, :]
#    d_a    = simdata.d_adjust[r0, :]
    adj    = simdata.adjust_indicator[r0, :]

 #   a_lag  = simdata.a[r1, :]
  #  aa_lag = simdata.aa[r1, :]
   # d_lag  = simdata.d[r1, :]
   # ex_lag = simdata.ex[r1, :]

    # --- value/denominators in local units, panel-wide ---
    Vd    = pd .* ex .* d                         # durable value
    a_loc = aa
    a_fx  = ex .* a                               # USD assets valued in local
    a_eff = a_loc .+ a_fx                         # total liquid in local
  #  income = max.(w .* y, 1e-12)                  # guard zeros

    # --- build panel vectors (flatten time×households) ---
    d_value_vec     = vec(Vd)
    a_eff_vec       = vec(a_eff)
    a_fx_vec        = vec(a_fx)
#    income_vec      = vec(income)

    # shares/ratios at (t,i)
   # d_income_ratio  = d_value_vec ./ income_vec
    d_wealth_ratio  = d_value_vec ./ max.(d_value_vec .+ a_eff_vec, 1e-12)
    usd_sh          = a_fx_vec    ./ max.(a_eff_vec, 1e-12)
    usd_pt          = usd_sh .> 0.0

    # --- durations at every (t,i) ---
    spells_periods = running_spells(adj)              # T×N
    duration_years = vec(spells_periods) ./ per_year  # flatten
   # adj_spell      = 1.0 ./ (1.0 .+ duration_years)

    # ---- moments over the entire panel ----
    m1 = mean(_fin(duration_years))
    m2 = mean(_fin(Float64.(usd_pt)))

    dwr_fin = _fin(d_wealth_ratio)
    m3 = mean(dwr_fin)
    m4 = var(dwr_fin)                     
    m5 = var(_fin(Float64.(usd_pt)))
 

    # Diagnostics you already use
  #  gap_vec, f_x, x_values, h_x, I_d_abs, mu_gap, var_gap, adj_rate_gap =
   #     adjustment_gaps_sim(d_lag, d_a, adj)

    outmoms = [m1, m2, m3, m4, m5]

    if any(.!isfinite.(outmoms))
        badix = findfirst(!isfinite, outmoms)
        @warn "makemoments (panel) produced NaN/Inf" bad_moment=badix bad_value=outmoms[badix]
    end

    if settings.verbose
        println("Moments over full panel (|r0|×N = $(length(r0))×$(size(adj,2)))")
        println("m1 duration = ", m1)
        println("m2 USD participation mean = ", m2)
        println("m3 d/wealth mean = ", m3)
        println("m4 d/wealth var  = ", m4)
        println("m5 USD participation var = ", m5)
    end

    return outmoms
end
