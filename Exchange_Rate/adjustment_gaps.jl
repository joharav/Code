using StatsBase
using KernelDensity
function adjustment_gaps(adjust_result::NamedTuple, noadjust_result::NamedTuple)
    # Extract grids and durable goods grid
    grids = adjust_result.g  # Assuming grids are part of the result
    d_grid = grids.d         # Durable goods grid

    # Conditional assignments based on the indicator_matrix
    d_a = adjust_result.pol.d

    # Ensure no zero or negative values for logarithm
    d_a = max.(d_a, 1e-10)
    d_grid = max.(d_grid, 1e-10)

    # Initialize gap array
    gap = zeros(sz.np, sz.ne, sz.na, sz.nd)

    # Compute the gap x = log(d*) - log(d_{-1})
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    gap[ip, ie, ia, id] = log(d_a[ip, ie, ia, id]) - log(d_grid[id])
                end
            end
        end
    end

    # Flatten the gap array to a vector
    gap_vec = vec(gap)

    # Compute the distribution of gaps f(x)
    kd       = kde(gap_vec)
    f_x      = kd.density
    x_values = collect(kd.x)  # Convert to Vector{Float64}

    # Compute the adjustment hazard h(x) as a vector
    indicator_matrix = adjust_result.v .> noadjust_result.v
    adjustment_indicator = indicator_matrix
    adjustment_vec = vec(adjustment_indicator)

    # Identify unique gap values for comparison
    #unique_gaps = unique(gap_vec, dims=1)  # Note: May result in removal of floating-point precision concerns

    # Adjustment probabilities
    h_x = zeros(length(x_values))

    for i in 1:length(x_values)
        # Using tolerance to approximate matching gaps
        matches = abs.(gap_vec .- x_values[i]) .< sz.toler

        if any(matches)  # Ensure some states matched this gap
            # Probability of adjustment for this gap value
            h_x[i] = mean(adjustment_vec[matches])
        else
            h_x[i] = 0.0  # Default value if no matches, shouldn't typically occur with tolerance
        end
    end

    
    # Compute the aggregate durable expenditures I_d using sum of products

    I_d = sum(x_values .* h_x .* f_x)

    return gap, f_x, x_values, h_x, I_d
end
