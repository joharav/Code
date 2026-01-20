# ==========================================================================
# 4D MODEL: Simulation
# State: (e, y, w, d)
# Tracks: dollar share s, derived aa and a for moment computation
# ==========================================================================

@inline clamp_to_grid(x, g::AbstractVector{<:Real}) = min(max(x, first(g)), last(g))

@inline function draw_next(cdfrow::AbstractVector{<:Real}, u::Real)
    @inbounds begin
        n = length(cdfrow)
        uu = min(max(u, 0.0), 1.0 - eps())
        j = searchsortedfirst(cdfrow, uu)
        return j < 1 ? 1 : (j > n ? n : j)
    end
end


function simmodel(answ::NamedTuple)
    # Unpack
    grids = answ.g
    tmat = grids.t
    exg = grids.ex
    yg = grids.y
    w_grid = grids.w
    d_grid = grids.d
    s_grid = grids.s
    
    # Parameters for wealth evolution
    beta = answ.pea[1]
    rr = (1 / beta) - 1
    rr_star = answ.pea[9]
    kappa = answ.pea[11]

    # Precompute joint (e,y) CDF
    phatcdf = cumsum(tmat, dims=2)
    @inbounds phatcdf[:, end] .= 1.0

    # Ergodic initialization
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0

    # Storage
    allv = zeros(sz.nYears, sz.nFirms)
    allw = zeros(sz.nYears, sz.nFirms)     # total wealth
    alla = zeros(sz.nYears, sz.nFirms)     # dollar assets
    allaa = zeros(sz.nYears, sz.nFirms)    # peso assets
    alls = zeros(sz.nYears, sz.nFirms)     # dollar share
    alle = zeros(sz.nYears, sz.nFirms)
    ally = zeros(sz.nYears, sz.nFirms)
    alld = zeros(sz.nYears, sz.nFirms)
    allc = zeros(sz.nYears, sz.nFirms)
    adjust_indicator = zeros(sz.nYears, sz.nFirms)

    # Initial states
    ls = zeros(Int, sz.nYears + 1, sz.nFirms)
    @inbounds for ifi in 1:sz.nFirms
        ls[1, ifi] = draw_next(cdf_wgt, globals.draws[1, ifi])
    end

    # Initial positions from random draws
    wstart = globals.draws[1, :]
    dstart = globals.draws[2, :]
    sstart = globals.draws[3, :]  # initial dollar share

    Threads.@threads for ifi in 1:sz.nFirms
        @inbounds begin
            # Initial state indices
            pickw = min(Int(floor(sz.nw * wstart[ifi])) + 1, sz.nw)
            pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
            
            pickey = ls[1, ifi]
            picke = div(pickey - 1, sz.ny) + 1
            picky = mod(pickey - 1, sz.ny) + 1
            
            # Initial levels
            w_old = w_grid[pickw]
            d_old = d_grid[pickd]
            s_old = sstart[ifi]  # initial dollar share
            
            # Derive initial aa and a
            E_old = exg[picke]
            aa_old = (1.0 - s_old) * w_old
            a_old = (s_old * w_old) / E_old

            for iti in 1:sz.nYears
                # Current aggregate states
                e = exg[picke]
                y = yg[picky]

                # Clamp continuous states
                w_in = clamp_to_grid(w_old, w_grid)
                d_in = clamp_to_grid(d_old, d_grid)

                # Compare adjustment vs non-adjustment values
                vA = interpol_ey(picke, picky, w_in, d_in, grids, answ.adjust_result.v)
                vN = interpol_ey(picke, picky, w_in, d_in, grids, answ.noadjust_result.v)
                do_adjust = (vA - vN) > sz.toler

                # Get policies from winning regime
                if do_adjust
                    w_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.adjust_result.pol.w)
                    d_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.adjust_result.pol.d)
                    s_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.adjust_result.pol.s)
                    c_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.adjust_result.pol.c)
                    v_pr = vA
                else
                    w_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.noadjust_result.pol.w)
                    d_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.noadjust_result.pol.d)
                    s_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.noadjust_result.pol.s)
                    c_pr = interpol_ey(picke, picky, w_in, d_in, grids, answ.noadjust_result.pol.c)
                    v_pr = vN
                end

                # Derive aa and a from w and s
                aa_pr = (1.0 - s_pr) * w_pr
                a_pr = (s_pr * w_pr) / max(e, 1e-10)

                # Record adjustment
                adjust_indicator[iti, ifi] = (do_adjust && abs(d_pr - d_in) > sz.toler) ? 1.0 : 0.0

                # Save paths
                allv[iti, ifi] = v_pr
                allw[iti, ifi] = w_pr
                alla[iti, ifi] = a_pr
                allaa[iti, ifi] = aa_pr
                alls[iti, ifi] = s_pr
                alle[iti, ifi] = e
                ally[iti, ifi] = y
                alld[iti, ifi] = d_pr
                allc[iti, ifi] = c_pr

                # Transition to next period
                # Draw next (e,y) state
                row = ls[iti, ifi]
                u = globals.draws[iti + 1, ifi]
                nxt = draw_next(view(phatcdf, row, :), u)
                ls[iti + 1, ifi] = nxt
                
                picke_new = div(nxt - 1, sz.ny) + 1
                picky_new = mod(nxt - 1, sz.ny) + 1
                e_new = exg[picke_new]

                # Wealth evolution with portfolio returns
                # w_{t+1} = (1-s)*w'*(1+r) + s*w'*(1+r*)*(E'/E) - κ*s*w'
                trans_cost = kappa * s_pr * w_pr
                w_new = (1.0 - s_pr) * w_pr * (1.0 + rr) + 
                       s_pr * w_pr * (1.0 + rr_star) * (e_new / e) - trans_cost
                w_new = max(w_new, w_grid[1])
                w_new = min(w_new, w_grid[end])

                # Update states
                picke = picke_new
                picky = picky_new
                w_old = w_new
                d_old = d_pr
                s_old = s_pr
                aa_old = aa_pr
                a_old = a_pr
            end
        end
    end

    return (
        v = allv,
        w = allw,
        d = alld,
        a = alla,
        aa = allaa,
        s = alls,
        ex = alle,
        y = ally,
        c = allc,
        adjust_indicator = adjust_indicator,
    )
end
