# ==========================================================================
# 4D MODEL: Bellman optimization for ADJUSTMENT regime
# 
# Key innovation: Portfolio choice s* is solved within-period for each (w', d')
# 
# Given w' saved for next period, household chooses s ∈ [0,1] (dollar share)
# This affects:
#   1. Transaction cost: κ * s * w'
#   2. Next period wealth: w_{t+1} = (1-s)*w'*(1+r) + s*w'*(1+r*)*(E'/E)
#
# The optimal s* depends on the exchange rate process and risk aversion.
# ==========================================================================

function maxbellman(queue::Array{Float64,4}, util::Array{Float64,6}, 
                   grids::NamedTuple, pea::Vector{Float64})
    
    beta = pea[1]
    rr = (1 / beta) - 1       # peso rate
    rr_star = pea[9]          # dollar rate
    kappa = pea[11]           # dollar transaction cost
    gamma = pea[6]            # risk aversion
    
    e_grid = grids.ex
    w_grid = grids.w
    wp_grid = grids.wp
    s_grid = grids.s
    trans_e = grids.te        # ER transition matrix (ne × ne)
    trans_joint = grids.t     # joint (e,y) transition
    
    # Output arrays
    vnew = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    gidx = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),  # w' index
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),  # d' index
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)   # s index (optimal portfolio)
    )

    Threads.@threads for id in 1:sz.nd
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    vstar = -1e10
                    wstar = 1
                    dstar = 1
                    sstar = 1
                    
                    E_now = e_grid[ie]
                    
                    # Search over wealth and durable policies
                    for idp in 1:sz.npd
                        for iwp in 1:sz.npw
                            # Flow utility from this (w', d') choice
                            u_flow = util[ie, iy, iw, id, iwp, idp]
                            
                            if u_flow > -1e9  # feasible choice
                                w_next = wp_grid[iwp]
                                
                                # Find optimal portfolio choice s* for this w'
                                best_s_idx = 1
                                best_EV = -Inf
                                
                                for is in 1:sz.ns
                                    s = s_grid[is]
                                    
                                    # Transaction cost (peso terms)
                                    trans_cost = kappa * s * w_next
                                    
                                    # Expected continuation value over ER states
                                    EV = 0.0
                                    
                                    for ie_next in 1:sz.ne
                                        E_next = e_grid[ie_next]
                                        
                                        # Wealth realization given portfolio
                                        # Peso part: (1-s)*w'*(1+r)
                                        # Dollar part: s*w'*(1+r*)*(E'/E)
                                        w_realized = (1.0 - s) * w_next * (1.0 + rr) + 
                                                    s * w_next * (1.0 + rr_star) * (E_next / E_now) -
                                                    trans_cost
                                        
                                        w_realized = max(w_realized, w_grid[1])
                                        w_realized = min(w_realized, w_grid[end])
                                        
                                        # Interpolate continuation value
                                        # queue is expected value over y', already computed
                                        # queue[ie', iy, iw, id'] where we need to interpolate iw
                                        
                                        iw_L, iw_U, wt_w = bracket_grid(w_realized, w_grid)
                                        
                                        # Get continuation values at bracket points
                                        # Note: queue already has E_y[V] built in from the 
                                        # trans_joint multiplication done before calling maxbellman
                                        V_L = queue[ie_next, iy, iw_L, idp]
                                        V_U = queue[ie_next, iy, iw_U, idp]
                                        V_interp = (1.0 - wt_w) * V_L + wt_w * V_U
                                        
                                        # Weight by ER transition probability
                                        prob_e = trans_e[ie, ie_next]
                                        EV += prob_e * V_interp
                                    end
                                    
                                    if EV > best_EV
                                        best_EV = EV
                                        best_s_idx = is
                                    end
                                end
                                
                                # Bellman equation
                                bellman = u_flow + beta * best_EV
                                
                                if bellman > vstar
                                    vstar = bellman
                                    wstar = iwp
                                    dstar = idp
                                    sstar = best_s_idx
                                end
                            end
                        end
                    end
                    
                    vnew[ie, iy, iw, id] = vstar
                    gidx.w[ie, iy, iw, id] = wstar
                    gidx.d[ie, iy, iw, id] = dstar
                    gidx.s[ie, iy, iw, id] = sstar
                end
            end
        end
    end
    
    return vnew, gidx
end


# ==========================================================================
# Helper: bracket a value in a sorted grid
# ==========================================================================
@inline function bracket_grid(x::Float64, grid::Vector{Float64})
    n = length(grid)
    @inbounds begin
        if x <= grid[1]
            return 1, 1, 0.0
        elseif x >= grid[n]
            return n, n, 0.0
        else
            j = searchsortedlast(grid, x)
            xL = grid[j]
            xU = grid[j+1]
            wt = (xU == xL) ? 0.0 : (x - xL) / (xU - xL)
            return j, j+1, wt
        end
    end
end
