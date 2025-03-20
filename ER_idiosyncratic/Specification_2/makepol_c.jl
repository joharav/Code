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
    rr          = (1 / beta) - 1         # Discount rate

    a = grid.a       # Asset grid
    d = grid.d       # Durable goods grid
    e = grid.ex       # Exchange rate grid
    y = grid.y       # Income grid



    cpol = zeros(sz.ne,sz.ny,sz.na,sz.nd);
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na;
            Threads.@threads for iy in 1:sz.ny;
                Threads.@threads for ie in 1:sz.ne;
                    if ind==0
                        cpol[ie,iy,ia,id] = w * h * (1-tau) * y[iy] +  a[ia] * (1 .+  rr) - e[ie] * pd * delta * chi * d[id] -   apol[ie,iy,ia,id]
                    else
                        cpol[ie,iy,ia,id] = w * h * (1-tau) * y[iy] +  a[ia] * (1 .+  rr) + e[ie] * pd * (1 - f) * (1 - delta) * d[id] -    apol[ie,iy,ia,id] - e[ie] * pd * dpol[ie,iy,ia,id] - w * h * ft * y[iy]
                    end;
                    cpol[ie,iy,ia,id] = max(cpol[ie,iy,ia,id], 0.1)  # Ensure c is at least 0.1

                end
            end
        end;
    end;
    return cpol::Array{Float64};
end;
