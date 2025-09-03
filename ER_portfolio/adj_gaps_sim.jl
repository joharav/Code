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
    
    # Define bins using the KDE x_values
    n_bins = length(x_values)  # Use KDE points as bin centers
    bin_edges = range(minimum(x_values), stop=maximum(x_values), length=n_bins+1)
    bin_indices = [searchsortedlast(bin_edges, g) for g in gap_vec]

    # Compute empirical hazard per bin
    hazard_empirical = zeros(n_bins)
    bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in 1:n_bins]

    for i in 1:n_bins
        in_bin = (bin_indices .== i)
        total_in_bin = sum(in_bin)
        adjusted_in_bin = sum(adjustment_indicator[in_bin])  # Sum since it's binary (1 = adjusted)

        if total_in_bin > 0
            hazard_empirical[i] = adjusted_in_bin / total_in_bin
        else
            hazard_empirical[i] = NaN  # Avoid division by zero
        end
    end

    num_adjustments = sum(adjustment_indicator.==1)
    total_states = length(adjustment_indicator)
    adjustment_ratio = num_adjustments / total_states
    
    println("Total Adjustments: $num_adjustments / $total_states")
    println("Adjustment Ratio: $(round(adjustment_ratio * 100, digits=2))%")

    # Compute h(x) = f_adj(x) / f(x), avoiding division by zero
    h_x = hazard_empirical

    # Compute the aggregate durable expenditures I_d
    I_d = sum(x_values .* h_x .* f_x)

    return gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio
end
