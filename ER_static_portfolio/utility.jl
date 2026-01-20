# ==========================================================================
# 4D MODEL: Utility for ADJUSTMENT regime
# State:  (e, y, w, d) where w is beginning-of-period liquid resources
# Choice: (w', d')  (portfolio choice s handled in continuation)
# util[ie,iy,iw,id,iwp,idp] = u(c, d')
# ==========================================================================

function utility(grids::NamedTuple, pea::Vector{Float64})
    w_grid, d_grid, wp_grid, dp_grid = grids.w, grids.d, grids.wp, grids.dp
    e_grid, y_grid = grids.ex, grids.y

    beta  = pea[1]
    delta = pea[2]
    nu    = pea[5]
    gamma = pea[6]

    f    = pea[7]
    wage = pea[8]
    pd   = pea[10]
    tau  = pea[12]
    h    = pea[13]
    ft   = pea[17]   # time cost scale

    util = fill(-1e10, sz.ne, sz.ny, sz.nw, sz.nd, sz.npw, sz.npd)

    Threads.@threads for idp in 1:sz.npd
        for iwp in 1:sz.npw, id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            E     = e_grid[ie]
            Y     = y_grid[iy]
            w_now = w_grid[iw]
            d_now = d_grid[id]

            w_next = wp_grid[iwp]
            d_next = dp_grid[idp]

            labor_income = Y * wage * h * (1.0 - tau)

            # If you allow adjustment: sell depreciated old durable (net of f)
            sale_value = E * pd * (1.0 - f) * (1.0 - delta) * d_now

            # Buy new durable
            durable_purchase = E * pd * d_next

            # Time/resource cost of adjusting (your current formula: wage*h*ft*Y)
            time_cost = wage * h * ft * Y

            c = labor_income + w_now + sale_value - durable_purchase - w_next - time_cost

            if c > 0
                d_eff = max(d_next, 1e-8)
                util[ie, iy, iw, id, iwp, idp] =
                    ((c^nu * d_eff^(1.0 - nu))^(1.0 - gamma)) / (1.0 - gamma)
            end
        end
    end

    if settings.verbose
        bad = count(==( -1e10), util)
        total = length(util)
        println("Utility (adjust): penalized = ", bad, " / ", total,
                " (", round(100 * bad / total, digits=1), "%)")
    end

    return util
end

