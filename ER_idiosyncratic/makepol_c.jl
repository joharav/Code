function makepol_c(apol::Array{Float64}, dpol::Array{Float64}, grid::NamedTuple, ind::Int64)
    beta  = pea[1]
    delta = pea[2]
    f     = pea[7]
    w     = pea[8]
    chi   = pea[9]
    pd    = pea[10]
    ft    = pea[11]
    tau   = pea[12]
    h     = pea[13]
    rr    = (1 / beta) - 1

    a = grid.a
    d = grid.d
    e = grid.ex
    y = grid.y

    cpol = zeros(sz.ne, sz.ny, sz.na, sz.nd)

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for iy in 1:sz.ny
                Threads.@threads for ie in 1:sz.ne
                    income = w * h * (1 - tau) * y[iy] + e[ie] * a[ia] * (1 + rr)

                    if ind == 0
                        # Non-adjust: keep durable, no purchase/sale/ft
                        c = income - e[ie] * apol[ie, iy, ia, id]
                    else
                        # Adjust: sell old, buy new, pay fixed time cost
                        sale_value = e[ie] * pd * (1 - f) * (1 - delta) * d[id]
                        durable_purchase = e[ie] * pd * dpol[ie, iy, ia, id]
                        time_cost = w * h * ft * y[iy]
                        c = income + sale_value - durable_purchase - e[ie] * apol[ie, iy, ia, id] - time_cost
                    end

                    cpol[ie, iy, ia, id] = max(c, 0.1)
                end
            end
        end
    end
    return cpol
end
