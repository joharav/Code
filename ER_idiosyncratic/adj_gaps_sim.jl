using StatsBase
using KernelDensity
using Interpolations
using Plots

function adjustment_gaps_sim(current_d,d_a,adjustment_indicator)
    # Calculate the adjustement gaps
    gaps    = log.(d_a) .- log.(current_d)
    gap_vec = vec(gaps)
    # Compute KDE for only the adjusted cases
    @views adjusted_gaps = gap_vec[adjustment_indicator .== 1]


    # Calculate the moments of the gap 
    mu_gap  = mean(gap_vec)
    var_gap = var(gap_vec)

    # Densities and prob of adjustment
    # Compute the distribution of gaps f(x)
    kd       = kde(gap_vec)
    f_x      = kd.density
    x_values = collect(kd.x)  



    if !isempty(adjusted_gaps)  # Ensure we have valid data
        kde_f_adj = kde(adjusted_gaps)
        
        # Interpolate densities to match the x-values from kd
        f_adj_interp = LinearInterpolation(kde_f_adj.x, kde_f_adj.density, extrapolation_bc=Line())
        f_adj_x = f_adj_interp.(x_values)
    else
        f_adj_x = zeros(length(x_values))  # No adjustments occurred
    end
    
    num_adjustments = sum(adjustment_indicator.==1)
    total_states = length(adjustment_indicator)
    adjustment_ratio = num_adjustments / total_states
    
    println("Total Adjustments: $num_adjustments / $total_states")
    println("Adjustment Ratio: $(round(adjustment_ratio * 100, digits=2))%")

    # Compute h(x) = f_adj(x) / f(x), avoiding division by zero
    h_x = max.(f_adj_x ./ (f_x .+ 1e-6), 0)

    # Compute the aggregate durable expenditures I_d
    I_d = sum(x_values .* h_x .* f_x)

    return gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio
end
