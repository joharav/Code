# ==========================================================================
# 4D MODEL: Utility for NO-ADJUSTMENT regime
# State: (e, y, w, d) where w = total liquid wealth
# Policy: (w') only - durables depreciate deterministically
# ==========================================================================

function utility_noadjust(grids::NamedTuple, pea::Vector{Float64})
    w_grid, d, wp, dp, e, y = grids.w, grids.d, grids.wp, grids.dp, grids.ex, grids.y

    beta, delta, nu, gamma = pea[1], pea[2], pea[5], pea[6]
    wage, pd, kappa, tau, h = pea[8], pea[10], pea[11], pea[12], pea[13]
    rr = (1 / beta) - 1
    rr_star = pea[9]
    chi = pea[16]   # maintenance effectiveness

    # In non-adjust regime, durables depreciate with maintenance
    # d' = (1 - δ*(1-χ)) * d
    d_next_vec = (1.0 .- delta .* (1.0 .- chi)) .* d
    
    # Map each durable state to nearest policy grid point
    idp_map = [argmin(abs.(dp .- d_next_vec[id])) for id in 1:sz.nd]

    # Flow utility array: util[ie, iy, iw, id, iwp]
    # Note: idp is determined by id, so we don't loop over it
    util = zeros(sz.ne, sz.ny, sz.nw, sz.nd, sz.npw)

    Threads.@threads for iwp in 1:sz.npw
        for id in 1:sz.nd
            for iw in 1:sz.nw
                for iy in 1:sz.ny
                    for ie in 1:sz.ne
                        # Current state
                        E = e[ie]
                        Y = y[iy]
                        w_now = w_grid[iw]
                        d_now = d[id]
                        
                        # Next period wealth choice
                        w_next = wp[iwp]
                        
                        # Durables depreciate with maintenance
                        d_keep = d_next_vec[id]
                        
                        # Maintenance cost
                        maintenance = E * pd * delta * chi * d_now
                        
                        # Resources
                        income = Y * wage * h * (1 - tau) + w_now
                        
                        # Consumption (no adjustment costs, no durable trade)
                        c = income - w_next - maintenance
                        
                        if c > 0 && d_keep > 0
                            util[ie, iy, iw, id, iwp] = 
                                ((c^nu * d_keep^(1 - nu))^(1 - gamma)) / (1 - gamma)
                        else
                            util[ie, iy, iw, id, iwp] = -1e10
                        end
                    end
                end
            end
        end
    end

    if settings.verbose
        bad = count(u -> u == -1e10, util)
        total = length(util)
        println("Utility (no-adjust): penalized = ", bad, " / ", total,
                " (", round(100*bad/total, digits=1), "%)")
    end
    
    return util, idp_map, d_next_vec
end
