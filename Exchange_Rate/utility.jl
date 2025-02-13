function utility(grids::NamedTuple, pea::Vector{Float64})
    # Extract grids
    a = grids.a       # Asset grid
    d = grids.d       # Durable goods grid
    ap = grids.ap     # Future asset grid
    dp = grids.dp     # Future durable goods grid
    pd = grids.p      # Durable price grid
    e = grids.ex       # Exchange rate grid

    # Model parameters
    beta = pea[1]     # Discount factor
    delta = pea[2]    # Depreciation rate for durables
    nu = pea[7]       # Share parameter for nondurables
    gamma = pea[8]    # Risk aversion
    f = pea[9]        # Adjustment cost
    w = pea[10]       # Wage rate

    rr = (1 / beta) - 1 
    # Initialize utility array
    util = zeros(sz.np, sz.ne, sz.na, sz.nd, sz.npa, sz.npd)

    Threads.@threads for iid in 1:sz.npd
        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                Threads.@threads for ia in 1:sz.na
                        Threads.@threads for ie in 1:sz.ne
                            Threads.@threads for ip in 1:sz.np
                            # Calculate consumption and durable goods stock
                            c = w + e[ie] * a[ia] * (1 + rr) + e[ie] * pd[ip] * (1 - f) * (1 - delta) * d[id] -  e[ie] *  ap[iia] - e[ie] * pd[ip] * dp[iid]

                            # Check feasibility of consumption and durable goods stock
                            if c > 0 && dp[iid] > 0
                                # Calculate utility
                                util[ip, ie, ia, id, iia, iid] = (((c^nu) * (dp[iid]^(1 - nu)))^(1 - gamma)) / (1 - gamma)
                            else
                                # Apply penalty for infeasible choices
                                util[ip, ie, ia, id, iia, iid] = -1e10
                            end
                        end
                    end
                end
            end
        end
    end

    return util
end
