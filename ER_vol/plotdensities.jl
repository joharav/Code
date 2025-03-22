function plotdensities(x_values::Vector{Float64}, f_x::Vector{Float64}, variable_name::String; shock::Bool = settings.irfsshock)
    filename_suffix = shock ? "_shock" : ""

    # Plot the density
    p1 = plot(x_values, f_x, label="Density", xlabel="$variable_name", ylabel="Density of $variable_name", legend=false, color=:blue, xlims=(0, 10))

    # Ensure the output directory exists
    output_dir = "Output/Ratios"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Save the plot with a specific filename
    savefig(p1, joinpath(output_dir, "Density_ratios_$variable_name$filename_suffix.png"))
end