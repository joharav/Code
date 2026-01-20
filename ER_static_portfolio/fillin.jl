# ==========================================================================
# 4D MODEL: Fill in value function from state grids to policy grids
# obj: (ne, ny, nw, nd) -> objlong: (ne, ny, npw, npd)
# ==========================================================================

function fillin(obj::Array{Float64,4}, g::NamedTuple)
    objlong = zeros(sz.ne, sz.ny, sz.npw, sz.npd)
    
    # Thread over (iwp, idp) combinations
    Threads.@threads for J in 1:(sz.npw * sz.npd)
        iwp = 1 + (J - 1) % sz.npw
        idp = 1 + (J - 1) ÷ sz.npw
        
        # w-weights (state w -> policy wp)
        w_down = Int(floor((sz.nw - 1.0) * (iwp - 1.0) / (sz.npw - 1)) + 1)
        w_up = (iwp == sz.npw) ? w_down : min(w_down + 1, sz.nw)
        den_w = g.w[w_up] - g.w[w_down]
        ww = (iwp == sz.npw || den_w == 0.0) ? 1.0 : (g.wp[iwp] - g.w[w_down]) / den_w
        
        # d-weights (state d -> policy dp)
        d_down = Int(floor((sz.nd - 1.0) * (idp - 1.0) / (sz.npd - 1)) + 1)
        d_up = (idp == sz.npd) ? d_down : min(d_down + 1, sz.nd)
        den_d = g.d[d_up] - g.d[d_down]
        wd = (idp == sz.npd || den_d == 0.0) ? 1.0 : (g.dp[idp] - g.d[d_down]) / den_d
        
        # 2D bilinear over (w, d) for each (ie, iy)
        for iy in 1:sz.ny
            for ie in 1:sz.ne
                v00 = obj[ie, iy, w_down, d_down]
                v01 = obj[ie, iy, w_down, d_up]
                v10 = obj[ie, iy, w_up, d_down]
                v11 = obj[ie, iy, w_up, d_up]
                
                v0 = (1 - ww) * v00 + ww * v10
                v1 = (1 - ww) * v01 + ww * v11
                
                objlong[ie, iy, iwp, idp] = (1 - wd) * v0 + wd * v1
            end
        end
    end
    
    return objlong
end
