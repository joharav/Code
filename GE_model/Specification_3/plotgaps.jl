default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plotgaps(x_values::Vector{Float64}, f_x::Vector{Float64}, h_x::Vector{Float64}, gap_vec::Vector{Float64}; shock::Bool = settings.irfsshock)
    filename_suffix = shock ? "_shock" : ""

    # Ensure the output directory exists
    output_dir = "Specification_3/Output/Gaps"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Plot the distribution f(x) on the left y-axis
    p = plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft, color=:blue)

    # Plot the adjustment probability h(x) on the right y-axis
    pp=twinx()
    plot!(pp, x_values, h_x, label="h(x)", ylabel="Hazard", legend=:topright, yaxis=:right, color=:red)

    # Save the plot
    savefig(p, joinpath(output_dir, "Gaps_distr_probability$filename_suffix.png"))

    # Plot the gap_vec as a histogram
    hist=histogram(gap_vec, bins=50, xlabel="Gap values", ylabel="Frequency", title="Histogram of Gap Values", legend=false)
    savefig(hist, joinpath(output_dir, "Gap_histogram$filename_suffix.png"))

    # Plot the distribution f(x) separately
    p1 = plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft, color=:blue)
    savefig(p1, joinpath(output_dir, "Gaps_distr$filename_suffix.png"))
    @save joinpath(output_dir, "Gaps_distr.jld2") p1

    # Plot the adjustment probability h(x) separately
    p2 = plot(x_values, h_x, label="h(x)", xlabel="Durable gap, log points", ylabel="Hazard", legend=:topright, color=:red)
    savefig(p2, joinpath(output_dir, "Gaps_probability$filename_suffix.png"))
    @save joinpath(output_dir, "Gaps_probability.jld2") p2
end
