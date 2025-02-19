using Plots

function plotgaps(f_x, x_values, h_x) 


    # Plot the distribution f(x) and the adjustment probability h(x)
    plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft)
    plot!(x_values, h_x, label="h(x)", ylabel="Hazard", legend=:topright, yaxis=:right)
    savefig("Output/Gaps/Gaps_distr_probability.png")

    plot(x_values, f_x, label="f(x)", xlabel="Durable gap, log points", ylabel="Density", legend=:topleft)
    savefig("Output/Gaps/Gaps_distr.png")

    plot(x_values, h_x, label="h(x)", xlabel="Durable gap, log points", ylabel="Hazard", legend=:topright)
    savefig("Output/Gaps/Gaps_probability.png")

end