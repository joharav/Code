# ==========================================================================
# 4D MODEL: GIRF simulation (single ER shock at T_shock)
# Paired innovations come from globals.draws (same as baseline simmodel).
# State evolves on (e,y,w,d). Policies read off grids via nearest_index.
#
# Adds:
#   d_adjust[t,i] = dp[ idpA ] where idpA is the adjust-regime target index
#                  evaluated at the (possibly shocked) decision index ie_dec.
#
# Options:
#   shock_mult: multiply current e by this at impact (if shock_to_max=false)
#   shock_to_max: override e at impact to maximum(exg)
#   shock_before_decisions: if true, agents condition on shocked e when choosing policies at impact
# ==========================================================================

@inline function draw_next_rowcdf(cdfrow::AbstractVector{<:Real}, u::Real)
    return searchsortedfirst(cdfrow, u)
end

@inline function nearest_index(grid::AbstractVector{<:Real}, x::Real)
    j = searchsortedlast(grid, x)
    if j < 1
        return 1
    elseif j >= length(grid)
        return length(grid)
    else
        return (abs(x - grid[j]) ≤ abs(grid[j+1] - x)) ? j : (j+1)
    end
end

function simmodel_girf(answ::NamedTuple, T_shock::Int;
                       shock_mult::Float64=1.2,
                       shock_to_max::Bool=false,
                       shock_before_decisions::Bool=true)

    g     = answ.g
    tmat  = g.t                 # (ne*ny)×(ne*ny), ROW-stochastic
    exg   = g.ex
    yg    = g.y
    wgrid = g.w                 # STATE grid
    wp    = g.wp                # POLICY grid
    dgrid = g.d                 # STATE grid
    dp    = g.dp                # POLICY grid
    sgrid = g.s

    pea     = answ.pea
    beta    = pea[1]
    delta   = pea[2]
    f       = pea[7]
    wage    = pea[8]
    rr_star = pea[9]
    pd      = pea[10]
    kappa   = pea[11]
    tau     = pea[12]
    h       = pea[13]
    ft      = pea[17]
    rr      = (1/beta) - 1

    # CDF rows for joint transition
    phatcdf = cumsum(tmat, dims=2)
    @inbounds phatcdf[:, end] .= 1.0

    # crude initial distribution (same as baseline)
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0

    T, N = sz.nYears, sz.nFirms

    # outputs (match baseline simmodel fields + d_adjust)
    allw = zeros(T, N)
    alld = zeros(T, N)
    alls = zeros(T, N)
    alle = zeros(T, N)
    ally = zeros(T, N)
    allc = zeros(T, N)
    adj  = zeros(T, N)
    d_adjust = zeros(T, N)

    d_next_vec = answ.noadjust_result.d_next_vec

    # initial joint state index
    ls = zeros(Int, T+1, N)
    @inbounds for i in 1:N
        ls[1, i] = searchsortedfirst(cdf_wgt, globals.draws[1, i])
    end

    # initial continuous picks
    wstart = globals.draws[1, :]
    dstart = globals.draws[2, :]

    emax = maximum(exg)

    Threads.@threads for i in 1:N
        # initial states from uniform draws
        iw0 = clamp(Int(floor(sz.nw * wstart[i])) + 1, 1, sz.nw)
        id0 = clamp(Int(floor(sz.nd * dstart[i])) + 1, 1, sz.nd)

        w_old = wgrid[iw0]
        d_old = dgrid[id0]

        pickey = ls[1, i]
        ie = div(pickey - 1, sz.ny) + 1
        iy = mod(pickey - 1, sz.ny) + 1

        for t in 1:T
            e_base = exg[ie]
            y      = yg[iy]

            # shock only changes the level used in this period
            e = (t == T_shock) ? (shock_to_max ? emax : (e_base * shock_mult)) : e_base

            # map continuous state -> indices for DP objects
            iw = nearest_index(wgrid, w_old)
            id = nearest_index(dgrid, d_old)

            # decision ER index at impact
            ie_dec = (t == T_shock && shock_before_decisions) ? nearest_index(exg, e) : ie

            # --- NEW: always store adjust-regime target durable under (possibly shocked) ie_dec ---
            idpA = answ.adjust_result.gidx.d[ie_dec, iy, iw, id]
            d_adjust[t, i] = dp[idpA]

            # realized regime (value-based indicator from your solver wrapper)
            do_adj = answ.adjustment_indicator[ie_dec, iy, iw, id]
            adj[t, i] = do_adj ? 1.0 : 0.0

            if do_adj
                iwp = answ.adjust_result.gidx.w[ie_dec, iy, iw, id]
                idp = answ.adjust_result.gidx.d[ie_dec, iy, iw, id]
                is  = answ.adjust_result.gidx.s[ie_dec, iy, iw, id]

                w_pr = wp[iwp]
                d_pr = dp[idp]
                s_pr = sgrid[is]

                labor_income     = y*wage*h*(1.0 - tau)
                sale_value       = e*pd*(1.0 - f)*(1.0 - delta)*dgrid[id]
                durable_purchase = e*pd*d_pr
                time_cost        = wage*h*ft*y
                c = labor_income + wgrid[iw] + sale_value - durable_purchase - w_pr - time_cost
            else
                iwp = answ.noadjust_result.gidx.w[ie_dec, iy, iw, id]
                is  = answ.noadjust_result.gidx.s[ie_dec, iy, iw, id]

                w_pr = wp[iwp]
                s_pr = sgrid[is]
                d_pr = d_next_vec[id]

                labor_income = y*wage*h*(1.0 - tau)
                c = labor_income + wgrid[iw] - w_pr
            end

            s_pr = clamp(s_pr, 0.0, 1.0)

            # record
            allw[t, i] = w_pr
            alld[t, i] = d_pr
            alls[t, i] = s_pr
            alle[t, i] = e
            ally[t, i] = y
            allc[t, i] = c

            # next (e,y)
            row = ls[t, i]
            u   = globals.draws[t+1, i]
            nxt = draw_next_rowcdf(view(phatcdf, row, :), u)
            ls[t+1, i] = nxt

            ie_new = div(nxt - 1, sz.ny) + 1
            iy_new = mod(nxt - 1, sz.ny) + 1
            e_new  = exg[ie_new]

            # wealth transition uses shocked e in denominator at impact
            trans_cost = kappa * s_pr * w_pr
            w_new = (1.0 - s_pr)*w_pr*(1.0 + rr) +
                    s_pr*w_pr*(1.0 + rr_star)*(e_new / max(e, 1e-12)) -
                    trans_cost

            # update continuous + Markov indices
            w_old = w_new
            d_old = d_pr
            ie    = ie_new
            iy    = iy_new
        end
    end

    return (
        w = allw, d = alld, d_adjust = d_adjust,
        s = alls, ex = alle, y = ally, c = allc,
        adjust_indicator = adj
    )
end
