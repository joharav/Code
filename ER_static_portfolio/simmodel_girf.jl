# ==========================================================================
# 4D MODEL: Generalized Impulse Response simulation
# Simulates with an exchange rate shock at time T_shock
# ==========================================================================

function simmodel_girf(answ::NamedTuple, T_shock::Int)
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
    
    # Joint (e,y) transitions
    phatcdf = cumsum(tmat, dims=2)
    @inbounds phatcdf[:, end] .= 1.0
    
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0
    
    # Pre-compute state indices
    ls = zeros(Int, nT + 1, nI)
    for i in 1:nI
        ls[1, i] = draw_next(cdf_wgt, globals.draws[1, i])
    end
    
    # Initial states
    wstart = globals.draws[1, :]
    dstart = globals.draws[2, :]
    sstart = globals.draws[3, :]
    
    # Shock value (max exchange rate)
    shock_val = maximum(exg)
    
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
                # Exchange rate (shocked at T_shock)
                e = (t == T_shock) ? shock_val : exg[picke]
                y = yg[picky]
                
                # Clamp states
                w_in = clamp(w_old, w_grid[1], w_grid[end])
                d_in = clamp(d_old, d_grid[1], d_grid[end])
                
                # Compare adjust vs non-adjust
                vA = interpol_ey(picke, picky, w_in, d_in, grids, answ.adjust_result.v)
                vN = interpol_ey(picke, picky, w_in, d_in, grids, answ.noadjust_result.v)
                do_adjust = (vA - vN) > sz.toler
                
                # Get policies
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
                
                # Derive aa and a
                aa_pr = (1.0 - s_pr) * w_pr
                a_pr = (s_pr * w_pr) / max(e, 1e-10)
                
                # Record
                adjust_indicator[t, i] = (do_adjust && abs(d_pr - d_in) > sz.toler) ? 1.0 : 0.0
                allv[t, i] = v_pr
                allw[t, i] = w_pr
                alla[t, i] = a_pr
                allaa[t, i] = aa_pr
                alls[t, i] = s_pr
                alle[t, i] = e
                ally[t, i] = y
                alld[t, i] = d_pr
                allc[t, i] = c_pr
                
                # Transition
                row = ls[t, i]
                u = globals.draws[t + 1, i]
                nxt = draw_next(view(phatcdf, row, :), u)
                ls[t + 1, i] = nxt
                
                picke_new = div(nxt - 1, sz.ny) + 1
                picky_new = mod(nxt - 1, sz.ny) + 1
                e_new = (t + 1 == T_shock) ? shock_val : exg[picke_new]
                
                # Wealth evolution
                trans_cost = kappa * s_pr * w_pr
                w_new = (1.0 - s_pr) * w_pr * (1.0 + rr) + 
                       s_pr * w_pr * (1.0 + rr_star) * (e_new / e) - trans_cost
                w_new = clamp(w_new, w_grid[1], w_grid[end])
                
                picke = picke_new
                picky = picky_new
                w_old = w_new
                d_old = d_pr
                s_old = s_pr
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
