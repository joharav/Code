using Plots
default(fontfamily = "Computer Modern")  # Looks like LaTeX

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

    # Size - Histogram
    adjust_size = abs.(vec(d[adjust_indicator .== 1]) .- vec(d_state[adjust_indicator .== 1]))
    histogram(adjust_size, bins=80, xlabel="Adjustment Size", ylabel="Frequency", color=:red,legend=false)
    savefig(joinpath(output_dir, "Adjustment_Size.png"))

end 

