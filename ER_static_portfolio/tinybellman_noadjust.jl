# ==========================================================================
# 4D MODEL: Local search Bellman for NO-ADJUSTMENT regime
# ==========================================================================

function tinybellman_noadjust(queue::Array{Float64,4}, util::Array{Float64,5},
                             idp_map::Vector{Int}, old_gidx::dtp.Ipol,
                             grids::NamedTuple, pea::Vector{Float64})
    
    beta = pea[1]
    rr = (1 / beta) - 1
    rr_star = pea[9]
    kappa = pea[11]
    
    e_grid = grids.ex
    w_grid = grids.w
    wp_grid = grids.wp
    s_grid = grids.s
    trans_e = grids.te
    
    vnew = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    gidx = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)
    )

    Threads.@threads for id in 1:sz.nd
        idp = idp_map[id]  # Fixed durable in no-adjust
        
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    old_iwp = old_gidx.w[ie, iy, iw, id]
                    old_is = old_gidx.s[ie, iy, iw, id]
                    
                    iwp_lo = max(1, old_iwp - sz.pad)
                    iwp_hi = min(sz.npw, old_iwp + sz.pad)
                    is_lo = max(1, old_is - 3)
                    is_hi = min(sz.ns, old_is + 3)
                    
                    vstar = -Inf
                    wstar, sstar = old_iwp, old_is
                    E_now = e_grid[ie]
                    
                    @inbounds for iwp in iwp_lo:iwp_hi
                        u_flow = util[ie, iy, iw, id, iwp]
                        
                        if u_flow > -1e9
                            w_next = wp_grid[iwp]
                            
                            for is in is_lo:is_hi
                                s = s_grid[is]
                                trans_cost = kappa * s * w_next
                                
                                EV = 0.0
                                for ie_next in 1:sz.ne
                                    E_next = e_grid[ie_next]
                                    w_realized = (1.0 - s) * w_next * (1.0 + rr) + 
                                                s * w_next * (1.0 + rr_star) * (E_next / E_now) -
                                                trans_cost
                                    w_realized = clamp(w_realized, w_grid[1], w_grid[end])
                                    
                                    iw_L, iw_U, wt = bracket_grid(w_realized, w_grid)
                                    V_interp = (1.0 - wt) * queue[ie_next, iy, iw_L, idp] +
                                               wt * queue[ie_next, iy, iw_U, idp]
                                    EV += trans_e[ie, ie_next] * V_interp
                                end
                                
                                bellman = u_flow + beta * EV
                                if bellman > vstar
                                    vstar = bellman
                                    wstar = iwp
                                    sstar = is
                                end
                            end
                        end
                    end
                    
                    vnew[ie, iy, iw, id] = vstar
                    gidx.w[ie, iy, iw, id] = wstar
                    gidx.d[ie, iy, iw, id] = idp
                    gidx.s[ie, iy, iw, id] = sstar
                end
            end
        end
    end
    
    return vnew, gidx
end
