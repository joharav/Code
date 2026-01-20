# ==========================================================================
# 4D MODEL: Bellman optimization for NO-ADJUSTMENT regime
# 
# In non-adjust regime:
# - Durables depreciate deterministically: d' = (1 - δ(1-χ)) * d
# - Only choice is (w', s) - how much to save and portfolio allocation
# ==========================================================================

function maxbellman_noadjust(queue::Array{Float64,4}, util::Array{Float64,5},
                            idp_map::Vector{Int}, grids::NamedTuple, pea::Vector{Float64})
    
    beta = pea[1]
    rr = (1 / beta) - 1       # peso rate
    rr_star = pea[9]          # dollar rate
    kappa = pea[11]           # dollar transaction cost
    
    e_grid = grids.ex
    w_grid = grids.w
    wp_grid = grids.wp
    s_grid = grids.s
    trans_e = grids.te        # ER transition matrix
    
    # Output arrays (same 4D structure)
    vnew = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    gidx = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),  # w' index
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),  # d' index (from idp_map)
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)   # s index
    )

    Threads.@threads for id in 1:sz.nd
        # In non-adjust, d' is determined by d
        idp = idp_map[id]
        
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    vstar = -Inf
                    wstar = 1
                    sstar = 1
                    
                    E_now = e_grid[ie]
                    
                    # Only search over wealth policy (d' is fixed)
                    for iwp in 1:sz.npw
                        # Flow utility
                        u_flow = util[ie, iy, iw, id, iwp]
                        
                        if u_flow > -1e9  # feasible
                            w_next = wp_grid[iwp]
                            
                            # Find optimal portfolio s*
                            best_s_idx = 1
                            best_EV = -Inf
                            
                            for is in 1:sz.ns
                                s = s_grid[is]
                                
                                # Transaction cost
                                trans_cost = kappa * s * w_next
                                
                                # Expected continuation value
                                EV = 0.0
                                
                                for ie_next in 1:sz.ne
                                    E_next = e_grid[ie_next]
                                    
                                    # Wealth realization
                                    w_realized = (1.0 - s) * w_next * (1.0 + rr) + 
                                                s * w_next * (1.0 + rr_star) * (E_next / E_now) -
                                                trans_cost
                                    
                                    w_realized = max(w_realized, w_grid[1])
                                    w_realized = min(w_realized, w_grid[end])
                                    
                                    # Interpolate
                                    iw_L, iw_U, wt_w = bracket_grid(w_realized, w_grid)
                                    
                                    V_L = queue[ie_next, iy, iw_L, idp]
                                    V_U = queue[ie_next, iy, iw_U, idp]
                                    V_interp = (1.0 - wt_w) * V_L + wt_w * V_U
                                    
                                    prob_e = trans_e[ie, ie_next]
                                    EV += prob_e * V_interp
                                end
                                
                                if EV > best_EV
                                    best_EV = EV
                                    best_s_idx = is
                                end
                            end
                            
                            bellman = u_flow + beta * best_EV
                            
                            if bellman > vstar
                                vstar = bellman
                                wstar = iwp
                                sstar = best_s_idx
                            end
                        end
                    end
                    
                    vnew[ie, iy, iw, id] = vstar
                    gidx.w[ie, iy, iw, id] = wstar
                    gidx.d[ie, iy, iw, id] = idp  # determined by state
                    gidx.s[ie, iy, iw, id] = sstar
                end
            end
        end
    end
    
    return vnew, gidx
end
