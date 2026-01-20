function utility_noadjust(grids::NamedTuple, pea::Vector{Float64})
    w_grid, d_grid, wp_grid, dp_grid = grids.w, grids.d, grids.wp, grids.dp
    e_grid, y_grid = grids.ex, grids.y

    beta  = pea[1]
    delta = pea[2]
    nu    = pea[5]
    gamma = pea[6]

    wage = pea[8]
    tau  = pea[12]
    h    = pea[13]
    chi  = pea[16]   # effectiveness in depreciation ONLY

    d_next_vec = (1.0 .- delta .* (1.0 .- chi)) .* d_grid
    idp_map = [argmin(abs.(dp_grid .- d_next_vec[id])) for id in 1:sz.nd]

    util = fill(-1e10, sz.ne, sz.ny, sz.nw, sz.nd, sz.npw)

    Threads.@threads for iwp in 1:sz.npw
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            Y     = y_grid[iy]
            w_now = w_grid[iw]
            w_next = wp_grid[iwp]

            labor_income = Y * wage * h * (1.0 - tau)
            c = labor_income + w_now - w_next

            if c > 0
                d_eff = max(d_next_vec[id], 1e-8)
                util[ie, iy, iw, id, iwp] =
                    ((c^nu * d_eff^(1.0 - nu))^(1.0 - gamma)) / (1.0 - gamma)
            end
        end
    end

    if settings.verbose
        bad = count(==( -1e10), util)
        total = length(util)
        println("Utility (no-adjust): penalized = ", bad, " / ", total,
                " (", round(100*bad/total, digits=1), "%)")
    end

    return util, idp_map, d_next_vec
end
