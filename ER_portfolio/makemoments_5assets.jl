using Statistics

# ---------- helpers ----------
_fin(v) = v[isfinite.(v)]
@inline _safe_den(x) = max(x, 1e-12)

# Running time-since-last-adjustment (in periods) for each (t,i) in r0
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
    delta = pea[2]  # depreciation rate

    # --- post-burn-in window ---
    r0 = (sz.burnin-2):sz.nYears
    r1 = (sz.burnin-3):(sz.nYears-1)   # currently unused, but fine

    # slice panel
    a      = simdata.a[r0, :]
    aa     = simdata.aa[r0, :]
    d      = simdata.d[r0, :]
    d_lag  = simdata.d[r1, :]
    ex     = simdata.ex[r0, :]
    y      = simdata.y[r0, :]
    adj    = simdata.adjust_indicator[r0, :]

    # if you have these in simdata, we’ll use them for volatility
    c       = hasproperty(simdata, :c)         ? simdata.c[r0, :]         : zeros(size(a))
    d_adj   = hasproperty(simdata, :d_adjust)  ? simdata.d_adjust[r0, :]  : zeros(size(d))

    # --- value/denominators in local units, panel-wide ---
    Vd    = pd .* ex .* d                         # durable stock value
    a_loc = aa
    a_fx  = ex .* a                               # USD assets valued in local
    a_eff = a_loc .+ a_fx                         # total liquid in local

    # --- build panel vectors (flatten time×households) ---
    d_value_vec = vec(Vd)
    a_eff_vec   = vec(a_eff)
    a_fx_vec    = vec(a_fx)
    c_vec       = vec(c)
    d_adj_vec   = vec(d_adj)

    # shares/ratios at (t,i)
    d_wealth_ratio  = d_value_vec ./ max.(d_value_vec .+ a_eff_vec, 1e-12)
    usd_sh          = a_fx_vec    ./ max.(a_eff_vec, 1e-12)
    usd_pt          = usd_sh .> 0.0

    # --- durations at every (t,i) ---
    spells_periods = running_spells(adj)              # T×N
    duration_years = vec(spells_periods) ./ per_year  # flatten

    # ---- core GMM moments over the entire panel ----
    # m1: mean duration (years)
    m1 = mean(_fin(duration_years))

    # d/wealth
    dwr_fin = _fin(d_wealth_ratio)
    m3 = mean(dwr_fin)          # dwealth_mean
    m4 = var(dwr_fin)           # dwealth_var

    # adj_rate
    adj_rate = mean(vec(adj) .> 0.0)
    m6 = adj_rate

    #usd assets
    m8 = mean(_fin(usd_sh))
    m9 = var(_fin(usd_sh))
    usd_share_mean = m8

    # ---- extra mechanisms: volatilities (NOT used in GMM) ----
    # consumption volatility: std(log(c + eps))
    eps = 1e-6
    cons_vol    = std(log.(_fin(c_vec .+ eps)))

    # durable spending: price × ex × d_adjust
    d_spend     = d - d_lag*(1-delta)
    d_spend_vec = vec(d_spend)
    d_spend_vol = std((_fin(d_spend_vec .+ eps)))

    # asset volatility: std(log(a_eff + eps))
    a_eff_vol   = std(log.(_fin(a_eff_vec .+ eps)))


    # outmoms = [duration_mean, dwealth_mean, dwealth_var, adj_rate, dollar mean, dollar vol]
    #outmoms = [m1, m3, m4, m6, m8, m9, cons_vol, d_spend_vol, a_eff_vol]
    outmoms = [duration_mean, dwealth_mean, dwealth_var, adj_rate, usd_share_mean ,m9, cons_vol, d_spend_vol, a_eff_vol]



    if any(.!isfinite.(outmoms))
        badix = findfirst(!isfinite, outmoms)
        @warn "makemoments (panel) produced NaN/Inf" bad_moment=badix bad_value=outmoms[badix]
    end

    if settings.verbose
        println("Moments over full panel (|r0|×N = $(length(r0))×$(size(adj,2)))")
        println("m1 duration_mean              = ", m1)
        println("m3 dwealth_mean               = ", m3)
        println("m4 dwealth_var                = ", m4)
        println("m6 adj_rate                   = ", m6)
        println("m8 dollar_share               = ", m8)
        println("m9 dollar_share_vol           = ", m9)
        println("cons_vol (log c)              = ", cons_vol)
        println("d_spend_vol (log pd*ex*d_adj) = ", d_spend_vol)
        println("a_eff_vol (log a_eff)         = ", a_eff_vol)
    end

    return outmoms
end
