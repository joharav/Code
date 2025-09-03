function makepol_c_twoasset(aa_pol::Array{Float64}, a_pol::Array{Float64}, dpol::Array{Float64},
    grid::NamedTuple, ind::Int64)
    beta  = pea[1];  delta = pea[2];  f = pea[7];  w = pea[8]
    chi   = pea[9];  pd    = pea[10]; ft = pea[11]; tau = pea[12]; h = pea[13]
    rr    = (1 / beta) - 1

    a   = grid.a     # foreign state grid
    aa  = grid.aa    # local   state grid
    d   = grid.d
    e   = grid.ex
    y   = grid.y

    cpol = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for iaa in 1:sz.na
                Threads.@threads for iy in 1:sz.ny
                    Threads.@threads for ie in 1:sz.ne
                        E  = e[ie]; Y = y[iy]
                        inc_assets = (1 + rr) * (aa[iaa] + E*a[ia])
                        inc_labor  = Y * w * h * (1 - tau)
                        income = inc_assets + inc_labor

                        next_pay = aa_pol[ie,iy,iaa,ia,id] + E * a_pol[ie,iy,iaa,ia,id]

                        if ind == 0
                            maint = E * pd * delta * chi * d[id]
                            c = income - next_pay - maint
                        else
                            sale_value = E * pd * (1 - f) * (1 - delta) * d[id]
                            purchase   = E * pd * dpol[ie,iy,iaa,ia,id]
                            time_cost  = Y * w * h * ft
                            c = income + sale_value - purchase - next_pay - time_cost
                        end

                        cpol[ie,iy,iaa,ia,id] = max(c, 0.1)
                    end
                end
            end
        end
    end
    return cpol
end
