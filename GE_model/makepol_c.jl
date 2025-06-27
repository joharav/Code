function makepol_c(apol::Array{Float64}, dpol::Array{Float64}, grid::NamedTuple, ind::Int64);
    beta        = pea[1]        # Discount factor
    delta       = pea[2]        # Depreciation rate for durables
    f           = pea[7]        # Adjustment cost
    w           = pea[8]        # Wage rate
    chi         = pea[9]        # Required maintenance
    pd          = pea[10]       # Price of durables 
    ft          = pea[11]       # fixed cost on wage rate
    tau         = pea[12]       # tax rate
    h           = pea[13]       # hours worked
    theta       = pea[16]       # Dollar share
    R_star      = pea[17]       # Dollar return
    R = (1 / beta)

    a = grid.a       # Asset grid
    d = grid.d       # Durable goods grid
    e = grid.ex       # Exchange rate grid
    zz = grid.zz       # Exchange rate grid



    cpol = zeros(sz.nz,sz.ne,sz.na,sz.nd);
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na;
            Threads.@threads for ie in 1:sz.ne;
                Threads.@threads for iz in 1:sz.nz;
                    y = w * h * (1 - tau) * zz[iz]
                    a_income = a[ia] * ((1 - theta) * R + theta * R_star * e[ie])
                    a_cost   = apol[iz,ie,ia,id] * ((1 - theta) + theta * e[ie])


                    if ind==0
                        cpol[iz,ie,ia,id] = y + a_income - e[ie] * pd * delta * chi * d[id] -  a_cost
                    else
                        cpol[iz,ie,ia,id] = y * (1-ft) + a_income + e[ie] * pd * (1 - f) * (1 - delta) * d[id] -   a_cost - e[ie] * pd * dpol[iz,ie,ia,id] 
                    end;
                    cpol[iz,ie,ia,id] = max(cpol[iz,ie,ia,id], 0.0001)  # Ensure c is at least 0.1
                end
            end
        end;
    end;
    return cpol::Array{Float64};
end;
