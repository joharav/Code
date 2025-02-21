using Plots

function plotgaps(x_values::Vector{Float64}, f_x::Vector{Float64}, h_x::Vector{Float64},gap_vec::Vector{Float64})
    # Plot the distribution f(x) on the left y-axis
    p = plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft, color=:blue)

    # Plot the adjustment probability h(x) on the right y-axis
    pp=twinx()
    plot!(pp, x_values, h_x, label="h(x)", ylabel="Hazard", legend=:topright, yaxis=:right, color=:red)

    # Save the plot
    savefig(p, "ER_vol/Output/Gaps/Gaps_distr_probability.png")

    # Plot the gap_vec as a histogram
    histogram(gap_vec, bins=50, xlabel="Gap values", ylabel="Frequency", title="Histogram of Gap Values", legend=false)
    savefig("ER_vol/Output/Gap_histogram.png")
    # Plot the distribution f(x) separately
    p1 = plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft, color=:blue)
    savefig(p1, "ER_vol/Output/Gaps/Gaps_distr.png")

    # Plot the adjustment probability h(x) separately
    p2 = plot(x_values, h_x, label="h(x)", xlabel="Durable gap, log points", ylabel="Hazard", legend=:topright, color=:red)
    savefig(p2, "ER_vol/Output/Gaps/Gaps_probability.png")

    # Plot the gap_vec as a histogram
    histogram(gap_vec, bins=50, xlabel="Gap values", ylabel="Frequency", title="Histogram of Gap Values", legend=false)
    savefig("ER_vol/Output/Gaps/Gap_histogram.png")

end