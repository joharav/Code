# ==========================================================================
# 4D MODEL: Howard policy improvement for NO-ADJUSTMENT regime
# Given fixed policy, iterate on value function only
# ==========================================================================

function howard_noadjust(queue::Array{Float64,4}, util::Array{Float64,5},
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
    
    Threads.@threads for id in 1:sz.nd
        idp = idp_map[id]  # Fixed durable index in no-adjust
        
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    # Use fixed policy
                    iwp = old_gidx.w[ie, iy, iw, id]
                    is = old_gidx.s[ie, iy, iw, id]
                    
                    # Flow utility
                    u_flow = util[ie, iy, iw, id, iwp]
                    
                    # Expected continuation with portfolio
                    w_next = wp_grid[iwp]
                    s = s_grid[is]
                    E_now = e_grid[ie]
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
                    
                    @inbounds vnew[ie, iy, iw, id] = u_flow + beta * EV
                end
            end
        end
    end
    
    # Keep durable policy pinned
    for id in 1:sz.nd
        old_gidx.d[:, :, :, id] .= idp_map[id]
    end
    
    return vnew, old_gidx
end
