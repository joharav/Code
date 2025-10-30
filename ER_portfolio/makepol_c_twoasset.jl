function makepol_c_twoasset(aa_pol::Array{Float64}, a_pol::Array{Float64}, dpol::Array{Float64},
    grid::NamedTuple, ind::Int64, pea::Vector{Float64})
    beta  = pea[1];  delta = pea[2];  f = pea[7];  w = pea[8]
    pd    = pea[10]; kappa = pea[11];  tau = pea[12]; h = pea[13]
    rr    = (1 / beta) - 1
    rr_foreign = pea[9]; 

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
                        inc_assets = (1 + rr) *aa[iaa] + (1+rr_foreign) * E*a[ia]
                        inc_labor  = Y * w * h * (1 - tau)
                        income = inc_assets + inc_labor
                        dollar_cost = kappa*(E * a_pol[ie,iy,iaa,ia,id])

                        next_pay = aa_pol[ie,iy,iaa,ia,id] + E * a_pol[ie,iy,iaa,ia,id]

                        if ind == 0
                            c = income - next_pay - dollar_cost
                        else
                            sale_value = E * pd * (1 - f) * (1 - delta) * d[id]
                            purchase   = E * pd * dpol[ie,iy,iaa,ia,id]
                            c = income + sale_value - purchase - next_pay - dollar_cost
                        end

                        cpol[ie,iy,iaa,ia,id] = max(c, 1e-6)
                    end
                end
            end
        end
    end
    return cpol
end
