using Statistics
using StatsBase

# ---------- helpers ----------
# Safe percent change: 100*(x - xlag)/xlag with guards against 0/neg/NaN
safe_change(x, xlag) = 100 .* (x .- xlag) ./ max.(abs.(xlag), 1e-12)

# Safe log with floor
@inline safe_log(x; eps=1e-12) = log.(max.(x, eps))

# Keep only finite entries from one or two aligned vectors
_fin(v) = v[isfinite.(v)]
_fin2(x, y) = begin
    m = isfinite.(x) .& isfinite.(y)
    (x[m], y[m])
end

# Decile index (1..10) based on x (unweighted); returns Vector{Int}
function decile_index(x::AbstractVector{<:Real})
    xf = _fin(x)
    qs = quantile(xf, 0:0.1:1)   # 11 cutpoints
    idx = similar(x, Int)
    @inbounds for i in eachindex(x)
        v = x[i]
        idx[i] = isfinite(v) ? clamp(searchsortedlast(qs, v), 1, 10) : 0
    end
    return idx
end

# Mean by decile with simple fallback if a decile is empty
function mean_by_decile(y::AbstractVector{<:Real}, dec::Vector{Int})
    out = fill(NaN, 10)
    @inbounds for d in 1:10
        m = dec .== d
        if any(m)
            out[d] = mean(_fin(y[m]))
        end
    end
    # fill NaNs with nearest neighbor to keep diffs/corrs defined
    for d in 1:10
        if !isfinite(out[d])
            l = findlast(isfinite, out[1:d])
            rpos = findfirst(isfinite, out[d:end])
            r = rpos === nothing ? d - 1 : d + rpos - 1
            out[d] = l === nothing && r < d ? 0.0 : (l === nothing ? out[r] : out[l])
        end
    end
    return out
end

# ---------- moments ----------
function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool=false)
    β  = pea[1]
    w  = pea[8]
    pd = pea[10]

    # aligned windows (quarterly); assumes globals: sz.burnin, sz.nYears
    t0 = (sz.burnin - 2):sz.nYears
    t1 = (sz.burnin - 3):(sz.nYears - 1)

    # pull series (T×N)
    a   = @view simdata.a[t0, :];     aa  = @view simdata.aa[t0, :]
    d   = @view simdata.d[t0, :];     ex  = @view simdata.ex[t0, :]
    y   = @view simdata.y[t0, :]
    d_a = @view simdata.d_adjust[t0, :]           # target d*
    adj = @view simdata.adjust_indicator[t0, :]   # 1 if adjusted at t
    d_lag = @view simdata.d[t1, :]               # lagged d for spells
    # core objects
    a_eff  = aa .+ ex .* a
    Vd     = pd .* ex .* d
    income = max.(w .* y, 1e-12)

    # ----- d-to-income distribution -----
    dinc   = vec(Vd) ./ vec(income)
    dinc_f = _fin(dinc)
    m1 = mean(dinc_f)
    m2 = var(dinc_f)
    qd = quantile(dinc_f, [0.25, 0.50, 0.75]); m3, m4, m5 = qd...

    # ----- dollarization -----
    a_loc  = aa
    a_fx   = ex .* a
    a_tot  = max.(a_loc .+ a_fx, 1e-12)
    usd_sh = vec(a_fx) ./ vec(a_tot)
    usd_pt = usd_sh .> 0.0
    m9     = mean(usd_pt)

    dec_idx = decile_index(vec(income))
    bydec   = mean_by_decile(usd_sh, dec_idx)
    m10     = bydec[9] - bydec[3]

    # ----- adjustment rate & corrs -----
    adj_flag = vec(adj) .> 0
    m6 = mean(adj_flag)

    xs, ys = _fin2(adj_flag, dinc);   m7  = cor(xs, ys)
    xs, ys = _fin2(adj_flag, usd_sh); m8  = cor(xs, ys)
    xs, ys = _fin2(usd_sh, dinc);     m11 = cor(xs, ys)
    xs, ys = _fin2(usd_sh, vec(a_eff)); m12 = cor(xs, ys)

    # ----- d-to-wealth distribution -----
    wealth = max.(Vd .+ a_eff, 1e-12)
    dw     = vec(Vd) ./ vec(wealth)
    dw_f   = _fin(dw)
    m13 = mean(dw_f)
    m14 = var(dw_f)
    qdw  = quantile(dw_f, [0.25, 0.50, 0.75]); m15, m16, m17 = qdw...

    # ----- gap–hazard diagnostics -----
    gap_vec, f_x, x_vals, h_x, I_d_abs, mu_gap, var_gap, adj_rate_gap =
        adjustment_gaps_sim(d_lag, d_a, adj)
    m18, m19, m20 = mu_gap, var_gap, I_d_abs
    # (adj_rate_gap should be ~ m6, but keep it as a cross-check if you like)

    # -----  duration spells of durable -----
    spells = spells_from_panel(adj)                # lengths in *quarters*
    spells_f = _fin(spells)
    m21 = mean(spells_f) /4                      # mean spell length (yy)
    qs  = quantile(spells_f, [0.50, 0.75])     # median & p75
    m22, m23 = qs...

    # -----  dispersion -----

    d_dispersion = var(log1p.(vec(Vd)))
    m24 = d_dispersion;


    outmoms = [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m21,m24]

    # out_named = (
    #     d_inc_mean      = m1,
    #     d_inc_var       = m2,
    #     d_inc_p25       = m3,
    #     d_inc_p50       = m4,
    #     d_inc_p75       = m5,
    #     adj_rate        = m6,
    #     corr_adj_dinc   = m7,
    #     corr_adj_usdsh  = m8,
    #     usd_particip    = m9,
    #     usdsh_p90m_p30  = m10,
    #     corr_usdsh_dinc = m11,
    #     corr_usdsh_aeff = m12,
    #     d_wealth_mean   = m13,
    #     d_wealth_var    = m14,
    #     d_wealth_p25    = m15,
    #     d_wealth_p50    = m16,
    #     d_wealth_p75    = m17,
    #     spell_mean_y    = m21,
    #     d_dispersion    = m24)

    return outmoms
end