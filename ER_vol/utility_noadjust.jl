function utility_noadjust(grids::NamedTuple, pea::Vector{Float64}) 
    # Extract grids
    a = grids.a       # Asset grid
    d = grids.d       # Durable goods grid
    ap = grids.ap     # Future asset grid
    e = grids.ex      # Exchange rate grid
    # Model parameters
    beta = pea[1]     # Discount factor
    delta = pea[2]    # Depreciation rate for durables
    nu = pea[5]       # Share parameter for nondurables
    gamma = pea[6]    # Risk aversion
    w = pea[8]        # Wage rate
    chi = pea[9]      # Required maintenance
    pd = pea[10]      # Price of durables 


    rr = (1 / beta) - 1 

    # Initialize utility array
    util = zeros(sz.ne, sz.na, sz.nd, sz.npa, sz.npd)
    ddp= zeros(sz.npd)

    Threads.@threads for iid in 1:sz.npd
        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                ddp[iid] = (1 - delta*(1-chi)) * d[id]  # Durable goods next period
                Threads.@threads for ia in 1:sz.na
                    Threads.@threads for ie in 1:sz.ne
                        # Calculate consumption and durable goods stock
                        c = w + e[ie] * a[ia] * (1 + rr) - e[ie] * pd * delta * chi * d[id] - e[ie] * ap[iia]

                        # Check feasibility of consumption and durable goods stock
                        if c > 0 && ddp[iid] > 0
                            # Calculate utility
                            util[ie, ia, id, iia,iid] = (((c^nu) * (ddp[iid]^(1 - nu)))^(1 - gamma)) / (1 - gamma)
                        else
                            util[ie, ia, id, iia,iid] = -1e10
                        end
                    end
                end
            end
        end
    end
    
    return util
end