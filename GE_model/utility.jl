function utility(grids::NamedTuple, pea::Vector{Float64})
    # Extract grids
    a = grids.a       # Asset grid
    d = grids.d       # Durable goods grid
    ap = grids.ap     # Future asset grid
    dp = grids.dp     # Future durable goods grid
    e = grids.ex       # Exchange rate grid
    zz = grids.zz       # Idiosyncratic Income grid

    # Model parameters
    beta    = pea[1]        # Discount factor
    delta   = pea[2]        # Depreciation rate for durables
    nu      = pea[5]        # Share parameter for nondurables
    gamma   = pea[6]        # Risk aversion
    f       = pea[7]        # Adjustment cost
    w       = pea[8]        # Wage rate
    pd      = pea[10]       # durable price
    ft      = pea[11]       # fixed cost on wage rate
    tau     = pea[12]       # tax rate
    h       = pea[13]       # hours worked
    theta       = pea[16]       # Dollar share
    R_star      = pea[17]       # Dollar return
    R = (1 / beta)

    # Initialize utility array
    util = zeros(sz.nz, sz.ne, sz.na, sz.nd, sz.npa,sz.npd)

    Threads.@threads for iid in 1:sz.npd
        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                Threads.@threads for ia in 1:sz.na
                    Threads.@threads for ie in 1:sz.ne
                        Threads.@threads for iz in 1:sz.nz
                            # Calculate consumption and durable goods stock
                            y = w * h * (1-tau)* zz[iz]
                            a_income = a[ia] * ((1 - theta) * R + theta * R_star * e[ie])
                            a_cost   = ap[iia] * ((1 - theta) + theta * e[ie])
                            
                            durable_return = e[ie] * pd * (1 - f) * (1 - delta) * d[id]
                            durable_cost   = e[ie] * pd * dp[iid]
                            
                            c = y * (1 - ft) + a_income + durable_return - a_cost - durable_cost

                            # Check feasibility of consumption and durable goods stock
                            if c > 0 && dp[iid] > 0
                                # Calculate utility
                                util[iz, ie, ia, id, iia, iid] = (((c^nu) * (ddp^(1 - nu)))^(1 - gamma)) / (1 - gamma)
                            else
                                util[iz, ie, ia, id, iia, iid] = -1e10

                            end
                        end
                    end
                end
            end
        end
    end


    return util
end
