using Plots

function d_adjust_time_size(simdata)
    
    d                   = simdata.d[sz.burnin-2:sz.nYears, :]
    d_state             = simdata.d[sz.burnin-3:sz.nYears-1, :]
    d_adjust            = simdata.d_adjust[sz.burnin-2:sz.nYears, :]
    adjust_indicator    = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]

    # Ensure the output directory exists
    output_dir="Output/Aggregates"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Timing - Histogram
    timing_of_adjustments = sum(adjust_indicator, dims=2)[:]  # Sum across firms (dims=2 → time)
    histogram(timing_of_adjustments, bins=30, xlabel="Time", ylabel="Frequency",legend=false)
    savefig(joinpath(output_dir, "Adjustment_Timing.png"))

    # Size - Histogram (filtering out tiny changes)
    adjust_size_raw = abs.(vec(d[adjust_indicator .== 1]) .- vec(d_state[adjust_indicator .== 1]))
    adjust_size = adjust_size_raw[adjust_size_raw .> 1e-2]  # e.g., filter out adjustments smaller than 0.01

    histogram(adjust_size, bins=80, xlabel="Adjustment Size", ylabel="Frequency", color=:red, legend=false)
    savefig(joinpath(output_dir, "Adjustment_Size.png"))


end 


function plot_simulated_d_and_a_by_state(simdata)
    d = simdata.d[sz.burnin-2:sz.nYears, :]
    a = simdata.a[sz.burnin-2:sz.nYears, :]
    e = simdata.ex[sz.burnin-2:sz.nYears, :]
    y = simdata.y[sz.burnin-2:sz.nYears, :]

    output_dir = "Output/Distr_State"
    isdir(output_dir) || mkpath(output_dir)

    unique_e = unique(vec(e))
    unique_y = unique(vec(y))

    for ei in unique_e, yi in unique_y
        # Cartesian indices
        idx = findall((e .== ei) .& (y .== yi))

        # Extract values directly from 2D arrays using Cartesian indices
        d_vals = [d[I] for I in idx]
        a_vals = [a[I] for I in idx]

        # Plot durable histogram
        histogram(d_vals, bins=50, xlabel="Durables", ylabel="Density", normalize=true,
                  title="Durables | e = $ei, y = $yi", legend=false)
        savefig(joinpath(output_dir, "d_dist_e$(ei)_y$(yi).png"))

        # Plot asset histogram
        histogram(a_vals, bins=50, xlabel="Assets", ylabel="Density", normalize=true,
                  title="Assets | e = $ei, y = $yi", legend=false)
        savefig(joinpath(output_dir, "a_dist_e$(ei)_y$(yi).png"))
    end
end

