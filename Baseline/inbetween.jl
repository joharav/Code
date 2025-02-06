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
    kwgt_lo = Int.(zeros(grid_size))
    kwgt_hi = Int.(zeros(grid_size))
    kwgt_fr = zeros(grid_size)

    # Debug: Print grid sizes
    println("grid_size: ", grid_size, " grid_policy_size: ", grid_policy_size)

    # Interpolation loop
    for ik = 1:grid_size
        kdown = Int(floor((grid_policy_size - 1.0) * (ik - 1.0) / (grid_size - 1)) + 1)
        kdown = min(kdown, grid_policy_size)  # Ensure kdown does not exceed bounds

        if ik == grid_size
            kup = kdown
            kfrac = 1.0
        else
            kup = min(kdown + 1, grid_policy_size)  # Ensure kup is within bounds
            if grid[kup] == grid[kdown]
                kfrac = 0.0
            else
                kfrac = (grid_policy[ik] - grid[kdown]) / (grid[kup] - grid[kdown])
            end
        end

        # Assign interpolation weights and indices
        kwgt_lo[ik] = kdown
        kwgt_hi[ik] = kup
        kwgt_fr[ik] = kfrac
    end

    # Return a named tuple containing the interpolation weights and endpoints
    k_wgt = (lo = kwgt_lo::Vector{Int64}, hi = kwgt_hi::Vector{Int64}, w = kwgt_fr::Vector{Float64})
    return k_wgt::NamedTuple{(:lo, :hi, :w)}
end
