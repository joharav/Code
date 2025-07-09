function utility_noadjust(grids::NamedTuple, pea::Vector{Float64}) 
    # Extract grids
    a = grids.a       # Asset grid
    d = grids.d       # Durable goods grid
    ap = grids.ap     # Future asset grid
    dp = grids.dp     # Future durable grid
    e = grids.ex      # Exchange rate grid
    y = grids.y       # Income grid
    # Model parameters
    beta        = pea[1]        # Discount factor
    delta       = pea[2]        # Depreciation rate for durables
    nu          = pea[5]        # Share parameter for nondurables
    gamma       = pea[6]        # Risk aversion
    w           = pea[8]        # Wage rate
    chi         = pea[9]        # Required maintenance
    pd          = pea[10]       # Price of durables 
    tau         = pea[12]       # tax rate
    h           = pea[13]       # hours worked

    rr = (1 / beta) - 1 

    penalty_count = 0  # Counter for penalties

    # Initialize utility array
    util = zeros(sz.ne, sz.ny, sz.na, sz.nd, sz.npa,sz.npd)
    ddp= (1 - delta * (1 - chi)) .* d  # Apply non-adjustment    
    iid = [argmin(abs.(dp .- ddp[id])) for id in 1:sz.nd]

    Threads.@threads for iid in 1:sz.npd
        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                Threads.@threads for ia in 1:sz.na
                    Threads.@threads for iy in 1:sz.ny
                        Threads.@threads for ie in 1:sz.ne
                            # Calculate consumption and durable goods stock
                            c = y[iy]* w * h * (1-tau) + e[ie] * a[ia] * (1 + rr) - e[ie] * pd * delta * chi * d[id] - e[ie] * ap[iia]
                            ddp = (1 - delta* (1 - chi)) * d[id]
                            # Check feasibility of consumption and durable goods stock
                            if c > 0 && ddp > 0
                                # Calculate utility
                                util[ie, iy, ia, id, iia,iid] = (((c^nu) * (ddp^(1 - nu)))^(1 - gamma)) / (1 - gamma)
                            else
                                util[ie, iy, ia, id, iia,iid] = -1e10
                                penalty_count += 1  # Increment counter for penalties

                            end
                        end
                    end
                end
            end
        end
    end
     penalty_share = penalty_count / (sz.ne * sz.ny * sz.na * sz.nd * sz.npa * sz.npd)
     println("Number of penalized states NA: ", penalty_count)
     println("Share of penalized states NA: ", penalty_share)




return util 





end
