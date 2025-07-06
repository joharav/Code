function utility_noadjust(grids::NamedTuple, pea::Vector{Float64}) 
    # Extract grids
    a = grids.a       # Asset grid
    d = grids.d       # Durable goods grid
    ap = grids.ap     # Future asset grid
    dp = grids.dp     # Future durable grid
    e = grids.ex      # Exchange rate grid
    zz = grids.zz       # Idiosyncratic Income grid

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
    theta       = pea[16]       # Dollar share
    R_star      = pea[17]       # Dollar return
    R = (1 / beta)


    # Initialize utility array
    util = zeros(sz.nz, sz.ne, sz.na, sz.nd, sz.npa,sz.npd)
    ddp_vec = (1 - delta * (1 - chi)) .* d
    iid = [argmin(abs.(dp .- ddp_vec[id])) for id in 1:sz.nd]

        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                Threads.@threads for ia in 1:sz.na
                    Threads.@threads for ie in 1:sz.ne
                        Threads.@threads for iz in 1:sz.nz
                            # Calculate consumption and durable goods stock
                            y = w * h * (1-tau)* zz[iz]
                            a_income = a[ia] * ((1 - theta) * R + theta * R_star * e[ie])
                            a_cost   = ap[iia] * ((1 - theta) + theta * e[ie])
                            c = y + a_income - e[ie] * pd * delta * chi * d[id] - a_cost
                            dnext = ddp_vec[id]  # This is a scalar
                            # Check feasibility of consumption and durable goods stock
                            if c > 0 && dnext > 0
                                # Calculate utility
                                util[iz, ie, ia, id, iia, iid[id]] = (((c^nu) * (dnext^(1 - nu)))^(1 - gamma)) / (1 - gamma)
                            else
                                util[iz, ie, ia, id, iia, iid[id]] = -1e10

                            end
                        end
                    end
                end
            end
        end
    
   
return util 



end
