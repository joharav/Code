using Plots

function plotgaps(x_values::Vector{Float64},
                  f_x::Vector{Float64},
                  h_x::Vector{Float64},
                  gap_vec::Vector{Float64};
                  shock::Bool = false,
                  outdir::String = "Output/Gaps")

    isdir(outdir) || mkpath(outdir)
    suffix = shock ? "_shock" : ""

    # Combined: density + hazard on right axis
    p = plot(x_values, f_x,
             xlabel="Durable gap (log points)",
             ylabel="Density f(x)",
             legend=false)
    p2 = twinx()
    plot!(p2, x_values, h_x,
          seriestype=:scatter,
          marker=:circle,
          ylabel="Hazard h(x)",
          legend=false)
    savefig(p, joinpath(outdir, "Gaps_distr_probability$(suffix).png"))

    # Histogram of gaps
    hist = histogram(gap_vec, bins=50,
                     xlabel="Gap values",
                     ylabel="Frequency",
                     title="Histogram of durable gaps",
                     legend=false)
    savefig(hist, joinpath(outdir, "Gap_histogram$(suffix).png"))

    # Separate density
    pD = plot(x_values, f_x,
              xlabel="Durable gap (log points)",
              ylabel="Density f(x)",
              legend=false)
    savefig(pD, joinpath(outdir, "Gaps_distr$(suffix).png"))

    # Separate hazard
    pH = plot(x_values, h_x,
              seriestype=:scatter,
              marker=:circle,
              xlabel="Durable gap (log points)",
              ylabel="Hazard h(x)",
              legend=false)
    savefig(pH, joinpath(outdir, "Gaps_probability$(suffix).png"))

    return nothing
end


function plotgaps_shock(x_values::Vector{Float64}, f_x::Vector{Float64}, h_x::Vector{Float64},
                        x_values_shock::Vector{Float64}, f_x_shock::Vector{Float64}, h_x_shock::Vector{Float64};
                        outdir::String = "Output/Gaps")

    isdir(outdir) || mkpath(outdir)

    p = plot(x_values, f_x,
             xlabel="Durable gap (log points)",
             ylabel="Density f(x)",
             label="No shock")
    plot!(p, x_values_shock, f_x_shock, label="Shock")
    savefig(p, joinpath(outdir, "Gaps_distr_comparison.png"))

    p2 = plot(x_values, h_x,
              xlabel="Durable gap (log points)",
              ylabel="Hazard h(x)",
              label="No shock")
    plot!(p2, x_values_shock, h_x_shock, label="Shock")
    savefig(p2, joinpath(outdir, "Gaps_probability_comparison.png"))

    return nothing
end


function plotdensities(x_values::Vector{Float64}, f_x::Vector{Float64}, variable_name::String;
                       shock::Bool = false,
                       outdir::String = "Output/Ratios",
                       xlims_override = nothing)

    isdir(outdir) || mkpath(outdir)
    suffix = shock ? "_shock" : ""

    p = plot(x_values, f_x,
             xlabel=variable_name,
             ylabel="Density",
             legend=false)
    if xlims_override !== nothing
        xlims!(p, xlims_override)
    end

    savefig(p, joinpath(outdir, "Density_$(variable_name)$(suffix).png"))
    return nothing
end
