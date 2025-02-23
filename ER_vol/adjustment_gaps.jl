using StatsBase
using KernelDensity
using Interpolations
using Plots

function adjustment_gaps(adjust_result::NamedTuple, noadjust_result::NamedTuple)
    # Extract grids and durable goods grid
    grids = adjust_result.g  
    d_grid = grids.d         # Durable goods grid

    # Conditional assignments based on the indicator_matrix
    d_a = adjust_result.pol.d

    # Ensure no zero or negative values for logarithm
    d_a = max.(d_a, 1e-10)
    d_grid = max.(d_grid, 1e-10)

    # Initialize gap array
    gap = zeros(sz.ne, sz.na, sz.nd)

    # Compute the gap x = log(d*) - log(d_{-1})
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                gap[ ie, ia, id] = log(d_a[ie, ia, id]) - log(d_grid[id])
            end
        end
    end

    # Flatten the gap array to a vector
    gap_vec = vec(gap)

    # Compute the distribution of gaps f(x)
    kd       = kde(gap_vec)
    f_x      = kd.density
    x_values = collect(kd.x)  

    # Compute the adjustment indicator
    adjustment_indicator = vec(adjust_result.v .> noadjust_result.v)

    # Compute KDE for only the adjusted cases
    adjusted_gaps = gap_vec[adjustment_indicator .== 1]

    if !isempty(adjusted_gaps)  # Ensure we have valid data
        kde_f_adj = kde(adjusted_gaps)
        
        # Interpolate densities to match the x-values from kd
        f_adj_interp = LinearInterpolation(kde_f_adj.x, kde_f_adj.density, extrapolation_bc=Line())
        f_adj_x = f_adj_interp.(x_values)
    else
        f_adj_x = zeros(length(x_values))  # No adjustments occurred
    end
    
    num_adjustments = sum(adjust_result.v .> noadjust_result.v)
    total_states = length(adjust_result.v)
    adjustment_ratio = num_adjustments / total_states
    
    println("Total Adjustments: ", num_adjustments)
    println("Total States: ", total_states)
    println("Adjustment Ratio: ", adjustment_ratio)

    histogram(vec(adjust_result.v .> noadjust_result.v), bins=2, xlabel="Adjustment (0 = No, 1 = Yes)", ylabel="Frequency", title="Histogram of Adjustments", legend=false)
    savefig("Output/Gap_adjustment_matrix.png")


    # Compute h(x) = f_adj(x) / f(x), avoiding division by zero
    h_x = f_adj_x ./ (f_x .+ 1e-10)

    # Compute the aggregate durable expenditures I_d
    I_d = sum(x_values .* h_x .* f_x)

    return gap, f_x, x_values, h_x, I_d
end
