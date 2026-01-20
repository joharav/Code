# ==========================================================================
# 4D MODEL: Generalized Impulse Response simulation (ER shock at T_shock)
# State: (e, y, w, d). Policies depend on (ie, iy, w, d).
# Shock is implemented by overriding the ER index at t == T_shock.
# ==========================================================================

@inline clamp_to_grid(x, g::AbstractVector{<:Real}) = min(max(x, first(g)), last(g))

function simmodel_girf(answ::NamedTuple, T_shock::Int)
    grids  = answ.g
    tmat   = grids.t
    exg    = grids.ex
    yg     = grids.y
    w_grid = grids.w
    d_grid = grids.d

    beta    = answ.pea[1]
    rr      = (1 / beta) - 1
    rr_star = answ.pea[9]
    kappa   = answ.pea[11]

    nT, nI = sz.nYears, sz.nFirms

    # Storage
    allv = zeros(nT, nI)
    allw = zeros(nT, nI)
    alla = zeros(nT, nI)
    allaa = zeros(nT, nI)
    alls = zeros(nT, nI)
    alle = zeros(nT, nI)
    ally = zeros(nT, nI)
    alld = zeros(nT, nI)
    allc = zeros(nT, nI)
    adjust_indicator = zeros(nT, nI)

    # Joint (e,y) transitions (rows = current combined state)
    phatcdf = cumsum(tmat, dims=2)
    @inbounds phatcdf[:, end] .= 1.0

    # Ergodic init over combined (e,y)
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0

    # Precompute shock ER index
    shock_val = maximum(exg)
    shock_ie = argmax(exg)  # index of shock_val

    # Initial combined (e,y) indices
    ls = zeros(Int, nT + 1, nI)
    @inbounds for i in 1:nI
        ls[1, i] = draw_next(cdf_wgt, globals.draws[1, i])
    end

    # Initial continuous states from draws
    wstart = globals.draws[1, :]
    dstart = globals.draws[2, :]
    sstart = globals.draws[3, :]

    Threads.@threads for i in 1:nI
        @inbounds begin
            pickw = min(Int(floor(sz.nw * wstart[i])) + 1, sz.nw)
            pickd = min(Int(floor(sz.nd * dstart[i])) + 1, sz.nd)

            pickey = ls[1, i]
            picke = div(pickey - 1, sz.ny) + 1
            picky = mod(pickey - 1, sz.ny) + 1

            w_old = w_grid[pickw]
            d_old = d_grid[pickd]
            s_old = sstart[i]

            for t in 1:nT
                # --------- current (e,y) indices for decision rules ----------
                ie_t = picke
                iy_t = picky

                # apply the shock by overriding the ER index AT THIS PERIOD
                if t == T_shock
                    ie_t = shock_ie
                end

                e = exg[ie_t]
                y = yg[iy_t]

                # Clamp continuous states
                w_in = clamp_to_grid(w_old, w_grid)
                d_in = clamp_to_grid(d_old, d_grid)

                # Compare regimes using SHOCKED INDEX if t == T_shock
                vA = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.adjust_result.v)
                vN = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.noadjust_result.v)
                do_adjust = (vA - vN) > sz.toler

                # Policies from winning regime (also using ie_t)
                if do_adjust
                    w_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.adjust_result.pol.w)
                    d_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.adjust_result.pol.d)
                    s_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.adjust_result.pol.s)
                    c_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.adjust_result.pol.c)
                    v_pr = vA
                else
                    w_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.noadjust_result.pol.w)
                    d_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.noadjust_result.pol.d)
                    s_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.noadjust_result.pol.s)
                    c_pr = interpol_ey(ie_t, iy_t, w_in, d_in, grids, answ.noadjust_result.pol.c)
                    v_pr = vN
                end

                # Derive aa and a at current e
                aa_pr = (1.0 - s_pr) * w_pr
                a_pr  = (s_pr * w_pr) / max(e, 1e-10)

                # Record
                adjust_indicator[t, i] = (do_adjust && abs(d_pr - d_in) > sz.toler) ? 1.0 : 0.0
                allv[t, i]  = v_pr
                allw[t, i]  = w_pr
                alla[t, i]  = a_pr
                allaa[t, i] = aa_pr
                alls[t, i]  = s_pr
                alle[t, i]  = e
                ally[t, i]  = y
                alld[t, i]  = d_pr
                allc[t, i]  = c_pr

                # --------- draw next (e,y) combined state ----------
                row = ls[t, i]
                u = globals.draws[t + 1, i]
                nxt = draw_next(view(phatcdf, row, :), u)
                ls[t + 1, i] = nxt

                picke_new = div(nxt - 1, sz.ny) + 1
                picky_new = mod(nxt - 1, sz.ny) + 1

                # next e index is the drawn one (no forced shock at t+1 unless you want it)
                ie_tp1 = picke_new
                e_new = exg[ie_tp1]

                # Wealth evolution uses e_new / e_current, where e_current is shocked if t==T_shock
                trans_cost = kappa * s_pr * w_pr
                w_new = (1.0 - s_pr) * w_pr * (1.0 + rr) +
                        s_pr * w_pr * (1.0 + rr_star) * (e_new / max(e, 1e-12)) -
                        trans_cost
                w_new = clamp_to_grid(w_new, w_grid)

                # Update states (note: picke itself follows the Markov draw; shock only overrides ie_t inside period)
                picke = picke_new
                picky = picky_new
                w_old = w_new
                d_old = d_pr
                s_old = s_pr
            end
        end
    end

    return (
        v = allv, w = allw, d = alld,
        a = alla, aa = allaa, s = alls,
        ex = alle, y = ally, c = allc,
        adjust_indicator = adjust_indicator,
    )
end

