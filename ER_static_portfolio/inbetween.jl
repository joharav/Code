function inbetween(grids::NamedTuple, islog::Bool, grid_type::Symbol)
    # Choose the appropriate grid based on grid_type
    if grid_type == :a
        grid = grids.a
        grid_policy = grids.ap
        grid_size = sz.na
        grid_policy_size = sz.npa
    elseif grid_type == :d
        grid = grids.d
        grid_policy = grids.dp
        grid_size = sz.nd
        grid_policy_size = sz.npd
    else
        error("Invalid grid type. Use :a or :d.")
    end

    # Optionally take the log of the grid if 'islog' is true
    if islog
        grid = log.(grid)
        grid_policy = log.(grid_policy)
    end

    # Initialize interpolation weights and indices
    kwgt_lo = Int.(zeros(grid_policy_size))
    kwgt_hi = Int.(zeros(grid_policy_size))
    kwgt_fr = zeros(grid_policy_size)

    # Interpolation loop
    Threads.@threads for ik = 1:grid_policy_size
        kdown = Int(floor((grid_size - 1.0) * (ik - 1.0) / (grid_policy_size - 1)) + 1)
        if ik == grid_policy_size
            kup = kdown
            kfrac = 1.0
        else
            kup = kdown + 1
            kfrac = (grid_policy[ik] - grid[kdown]) / (grid[kup] - grid[kdown])
        end
        kwgt_lo[ik] = kdown
        kwgt_hi[ik] = kup
        kwgt_fr[ik] = kfrac
    end

    k_wgt = (lo = kwgt_lo, hi = kwgt_hi, w = kwgt_fr)
    return k_wgt
end
