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

function simmodel(answ::NamedTuple)
    g     = answ.g
    tmat  = g.t
    exg, yg = g.ex, g.y
    wgrid, wp = g.w, g.wp
    dgrid, dp = g.d, g.dp
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

    phatcdf = cumsum(tmat, dims=2)
    @inbounds phatcdf[:, end] .= 1.0

    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0

    T, N = sz.nYears, sz.nFirms
    allw = zeros(T, N); alld = zeros(T, N); alls = zeros(T, N)
    alle = zeros(T, N); ally = zeros(T, N); allc = zeros(T, N)

    move_ind = zeros(T, N)     # tenure-reset event
    d_adjust = zeros(T, N)

    d_next_vec = answ.noadjust_result.d_next_vec

    ls = zeros(Int, T+1, N)
    @inbounds for i in 1:N
        ls[1, i] = searchsortedfirst(cdf_wgt, globals.draws[1, i])
    end

    wstart = globals.draws[1, :]
    dstart = globals.draws[2, :]

    tol = 1e-6

    Threads.@threads for i in 1:N
        iw0 = clamp(Int(floor(sz.nw * wstart[i])) + 1, 1, sz.nw)
        id0 = clamp(Int(floor(sz.nd * dstart[i])) + 1, 1, sz.nd)

        w_old = wgrid[iw0]
        d_old = dgrid[id0]

        pickey = ls[1, i]
        ie = div(pickey - 1, sz.ny) + 1
        iy = mod(pickey - 1, sz.ny) + 1

        for t in 1:T
            e = exg[ie]
            y = yg[iy]

            iw = nearest_index(wgrid, w_old)
            id = nearest_index(dgrid, d_old)

            # counterfactual adjust target every period
            idpA = answ.adjust_result.gidx.d[ie, iy, iw, id]
            d_adjust[t, i] = dp[idpA]

            # no-adjust implied next durable (passive depreciation)
            d_keep = d_next_vec[id]

            do_adj = answ.adjustment_indicator[ie, iy, iw, id]

            if do_adj
                iwp = answ.adjust_result.gidx.w[ie, iy, iw, id]
                idp = answ.adjust_result.gidx.d[ie, iy, iw, id]
                is  = answ.adjust_result.gidx.s[ie, iy, iw, id]

                w_pr = wp[iwp]
                d_pr = dp[idp]
                s_pr = sgrid[is]

                labor_income     = y*wage*h*(1.0 - tau)
                sale_value       = e*pd*(1.0 - f)*(1.0 - delta)*dgrid[id]
                durable_purchase = e*pd*d_pr
                time_cost        = wage*h*ft*y
                c = labor_income + wgrid[iw] + sale_value - durable_purchase - w_pr - time_cost
            else
                iwp = answ.noadjust_result.gidx.w[ie, iy, iw, id]
                is  = answ.noadjust_result.gidx.s[ie, iy, iw, id]

                w_pr = wp[iwp]
                s_pr = sgrid[is]
                d_pr = d_keep

                labor_income = y*wage*h*(1.0 - tau)
                c = labor_income + wgrid[iw] - w_pr
            end

            # tenure reset = “changed house” event
            move = do_adj && (abs(d_pr - d_keep) > tol * max(1.0, abs(d_keep)))
            move_ind[t, i] = move ? 1.0 : 0.0

            s_pr = clamp(s_pr, 0.0, 1.0)

            allw[t, i] = w_pr
            alld[t, i] = d_pr
            alls[t, i] = s_pr
            alle[t, i] = e
            ally[t, i] = y
            allc[t, i] = c

            row = ls[t, i]
            u   = globals.draws[t+1, i]
            nxt = draw_next_rowcdf(view(phatcdf, row, :), u)
            ls[t+1, i] = nxt

            ie_new = div(nxt - 1, sz.ny) + 1
            iy_new = mod(nxt - 1, sz.ny) + 1
            e_new  = exg[ie_new]

            trans_cost = kappa * s_pr * w_pr
            w_new = (1.0 - s_pr)*w_pr*(1.0 + rr) +
                    s_pr*w_pr*(1.0 + rr_star)*(e_new / e) -
                    trans_cost

            w_old = w_new
            d_old = d_pr
            ie = ie_new
            iy = iy_new
        end
    end

    return (
        w = allw, d = alld, d_adjust = d_adjust,
        s = alls, ex = alle, y = ally, c = allc,
        adjust_indicator = move_ind
    )
end

