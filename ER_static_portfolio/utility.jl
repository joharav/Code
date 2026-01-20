# ==========================================================================
# 4D MODEL: Utility for ADJUSTMENT regime
# State: (e, y, w, d) where w = total liquid wealth
# Policy: (w', d') with within-period portfolio choice s*
# ==========================================================================

function utility(grids::NamedTuple, pea::Vector{Float64})
    w_grid, d, wp, dp, e, y = grids.w, grids.d, grids.wp, grids.dp, grids.ex, grids.y
  
    beta, delta, nu, gamma = pea[1], pea[2], pea[5], pea[6]
    f, wage, pd, kappa, tau, h = pea[7], pea[8], pea[10], pea[11], pea[12], pea[13]
    rr = (1 / beta) - 1      # peso interest rate
    rr_star = pea[9]         # dollar interest rate
    ft = pea[17]             # time cost parameter

    # Flow utility array: util[ie, iy, iw, id, iwp, idp]
    # This is flow utility BEFORE portfolio choice adjustment
    # Portfolio choice affects continuation value, computed separately in Bellman
    util = zeros(sz.ne, sz.ny, sz.nw, sz.nd, sz.npw, sz.npd)

    Threads.@threads for idp in 1:sz.npd
        for iwp in 1:sz.npw
            for id in 1:sz.nd
                for iw in 1:sz.nw
                    for iy in 1:sz.ny
                        for ie in 1:sz.ne
                            # Current state
                            E = e[ie]
                            Y = y[iy]
                            w_now = w_grid[iw]
                            d_now = d[id]
                            
                            # Next period choices
                            w_next = wp[iwp]
                            d_next = dp[idp]
                            
                            # Resources available:
                            # - Income: Y * wage * h * (1 - tau)
                            # - Wealth carrying forward depends on last period's 
                            #   portfolio choice s_{t-1}, which determined how w_now 
                            #   was split between peso and dollar assets
                            # 
                            # SIMPLIFICATION: In steady state / simulation, we track
                            # the composition. For utility computation, we assume
                            # wealth earns a weighted average return.
                            # The portfolio choice s* for NEXT period is what matters
                            # for continuation value, not for today's flow utility.
                            
                            # For flow utility, assume w_now already reflects 
                            # returns from previous period's allocation
                            # (this is handled in simulation via tracking s)
                            
                            income = Y * wage * h * (1 - tau) + w_now
                            
                            # Durable adjustment: sell old, buy new
                            # Sale value of depreciated durables
                            sale_value = E * pd * (1 - f) * (1 - delta) * d_now
                            
                            # Purchase price of new durables
                            durable_purchase = E * pd * d_next
                            
                            # Time cost of adjustment
                            time_cost = wage * h * ft * Y
                            
                            # Consumption (note: kappa cost handled in continuation value)
                            c = income + sale_value - durable_purchase - w_next - time_cost
                            
                            if c > 0 && d_next > 0
                                util[ie, iy, iw, id, iwp, idp] = 
                                    ((c^nu * d_next^(1 - nu))^(1 - gamma)) / (1 - gamma)
                            else
                                util[ie, iy, iw, id, iwp, idp] = -1e10
                            end
                        end
                    end
                end
            end
        end
    end

    if settings.verbose
        bad = count(u -> u == -1e10, util)
        total = length(util)
        println("Utility (adjust): penalized = ", bad, " / ", total, 
                " (", round(100*bad/total, digits=1), "%)")
    end
    
    return util
end


# ==========================================================================
# Optimal portfolio choice function
# Given w' and current ER state, find s* that maximizes expected value
# ==========================================================================

function optimal_portfolio(w_next::Float64, ie::Int, V_next::Array{Float64,4}, 
                          grids::NamedTuple, pea::Vector{Float64})
    # Parameters
    beta = pea[1]
    rr = (1 / beta) - 1      # peso rate
    rr_star = pea[9]         # dollar rate  
    kappa = pea[11]          # dollar transaction cost
    
    e_grid = grids.ex
    trans_e = grids.te       # ER transition matrix
    E_now = e_grid[ie]
    
    # Grid search over dollar share s ∈ [0, 1]
    s_grid = grids.s
    ns = length(s_grid)
    
    best_s = 0.0
    best_val = -Inf
    
    # For small w', don't bother with dollars (transaction cost dominates)
    if w_next < 1e-6
        return 0.0, 0.0
    end
    
    @inbounds for is in 1:ns
        s = s_grid[is]
        
        # Transaction cost for dollar holdings
        # s * w_next is the peso value going into dollars
        trans_cost = kappa * s * w_next
        
        # Expected continuation value
        E_V = 0.0
        
        for ie_next in 1:sz.ne
            E_next = e_grid[ie_next]
            
            # Next period's wealth given portfolio choice:
            # Peso portion: (1-s) * w_next grows at (1+rr)
            # Dollar portion: s * w_next / E_now dollars, 
            #                 grows to s * w_next / E_now * (1+rr_star) dollars,
            #                 worth s * w_next * (1+rr_star) * E_next / E_now pesos
            
            w_next_realized = (1 - s) * w_next * (1 + rr) + 
                              s * w_next * (1 + rr_star) * (E_next / E_now) - 
                              trans_cost
            
            # Clamp to grid
            w_next_realized = max(w_next_realized, grids.w[1])
            w_next_realized = min(w_next_realized, grids.w[end])
            
            # Interpolate V_next at this wealth level
            # V_next is indexed [ie, iy, iw, id] - we need to integrate over iy too
            # For portfolio choice, we're computing E[V | e'] so we need V(e', y', w', d')
            # But y' is uncertain too. The full expectation is done in the Bellman.
            # Here we just compute the e' part for a given (iy, id) which will be
            # summed over in the main Bellman loop.
            
            # Get bracket for interpolation
            iw_L, iw_U, wt = bracket_wealth(w_next_realized, grids.w)
            
            # Probability of this ER transition
            prob_e = trans_e[ie, ie_next]
            
            # We return the expected V contribution from this ER state
            # The caller will handle the y transition and sum over id
            E_V += prob_e  # placeholder - actual interpolation done in Bellman
        end
        
        # For now, simplified: just track best s based on expected return
        # Full implementation integrates over the joint (e',y') transition
        expected_return = (1 - s) * (1 + rr) + s * (1 + rr_star) * 1.0 - kappa * s
        
        if expected_return > best_val
            best_val = expected_return
            best_s = s
        end
    end
    
    return best_s, best_val
end


# Helper: bracket a value in a grid
@inline function bracket_wealth(w::Float64, w_grid::Vector{Float64})
    n = length(w_grid)
    if w <= w_grid[1]
        return 1, 1, 0.0
    elseif w >= w_grid[n]
        return n, n, 0.0
    else
        j = searchsortedlast(w_grid, w)
        wL = w_grid[j]
        wU = w_grid[j+1]
        wt = (wU == wL) ? 0.0 : (w - wL) / (wU - wL)
        return j, j+1, wt
    end
end
