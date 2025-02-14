function adjustment_gaps(pea::Vector{Float64})
    # First, we are gonna recover all the durable policies, for the case where we adjust and when we do not

    # Adjusted and non-adjusted value functions
    adjust_result = valfun_adjust(pea)
    noadjust_result = valfun_noadjust(pea)
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
    kde = kde(gap_vec)
    f_x = kde.density
    x_values = kde.x

    # Compute the adjustment hazard h(x)
    indicator_matrix = adjust_result.v .> noadjust_result.v
    adjustment_indicator = indicator_matrix
    adjustment_vec = vec(adjustment_indicator)
    h_x = mean(adjustment_vec)
    
    # Compute the aggregate durable expenditures I_d using sum of products
    I_d = sum(x_values .* h_x .* f_x)

    return gap, f_x, x_values, h_x, I_d
end
