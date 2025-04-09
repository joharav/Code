default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plotgaps_comp(x_values::Vector{Float64}, f_x::Vector{Float64}, h_x::Vector{Float64}, gap_vec::Vector{Float64}, param_name::String, param_value::Float64)
    # Create directory if it doesn't exist
    output_dir = "Specification_2/Output/Gaps/Comparative/"
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # Generate unique filenames based on parameter and value
    param_str = string(param_name, "_", @sprintf("%.4f", param_value))

    # Plot the distribution f(x) on the left y-axis
    p = plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft, color=:blue)

    # Plot the adjustment probability h(x) on the right y-axis
    pp = twinx()
    plot!(pp, x_values, h_x, label="h(x)", ylabel="Hazard", legend=:topright, yaxis=:right, color=:red)

    # Save the combined plot
    savefig(p, output_dir * "Gaps_distr_probability_$(param_str).png")

    # Plot the gap_vec as a histogram
    histogram(gap_vec, bins=50, xlabel="Gap values", ylabel="Frequency", title="Histogram of Gap Values", legend=false)
    savefig(output_dir * "Gap_histogram_$(param_str).png")

    # Save individual plots separately
    p11 = plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft, color=:blue)
    p22 = plot(x_values, h_x, label="h(x)", xlabel="Durable gap, log points", ylabel="Hazard", legend=:topright, color=:red)

    @load joinpath(output_dir, "Gaps_distr.jld2") p1
    @load joinpath(output_dir, "Gaps_probability.jld2") p2

    # Combine the plots
    p_combined = plot(p1)
    plot!(p_combined, p11)
    savefig(p_combined, joinpath(output_dir, "Combined_Distr_$(param_str).png"))

    p_combined2 = plot(p22)
    plot!(p_combined2, p22)
    savefig(p_combined2, joinpath(output_dir, "Combined_Prob_$(param_str).png"))

    println("Saved plots for $(param_name) = $(param_value) in Output/Gaps/")
end