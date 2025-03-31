using Plots

function plotgaps_shock(x_values::Vector{Float64}, f_x::Vector{Float64}, h_x::Vector{Float64},  
                  x_values_shock::Vector{Float64}, f_x_shock::Vector{Float64}, h_x_shock::Vector{Float64})
    
    # Ensure the output directory exists
    output_dir = "Output/Gaps"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Plot f(x) and f(x) with shock on the same graph
    p = plot(x_values, f_x, label="No Shock", xlabel="Durable gap, log points, x", ylabel="Density", legend=:topright, color=:blue)
    plot!(p, x_values_shock, f_x_shock, label="Shock", color=:darkred)
    savefig(p, joinpath(output_dir, "Gaps_distr_comparison.png"))

    # Plot h(x) and h(x) with shock on the same graph
    p2 = plot(x_values, h_x, label="h(x) No Shock", xlabel="Durable gap, log points, x", ylabel="Hazard", legend=:topright, color=:green)
    plot!(p2, x_values_shock, h_x_shock, label="h(x) Shock", color=:purple)
    savefig(p2, joinpath(output_dir, "Gaps_probability_comparison.png"))

end