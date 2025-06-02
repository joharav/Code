using Plots
default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plotstuff(vee::Array{Float64, 4}, apol::Array{Int64, 4}, dpol::Array{Int64, 4},  cpol::Array{Float64, 4}, g::NamedTuple)

    a = g.a  # Assets
    d = g.d  # Durables

    # Iterate over price and exchange rate slices (fix: a, d ordering)
    for ie in 1:sz.ne
        dvee_slice = vee[:, ie, :, :]  
        dapol_slice = apol[:, ie, :, :]
        ddpol_slice = dpol[:, ie, :, :]
        dcpol_slice = cpol[:, ie, :, :]

        plot1 = surface(a, d, dvee_slice, xlabel="Assets", ylabel="Durables", zlabel="Value Function",legend=false)
        savefig(plot1, "Output/Policy/vf_slice_e$(ie).png")

        plot2 = surface(a, d, dapol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets",legend=false)
        savefig(plot2, "Output/Policy/Apolicy_slice_e$(ie).png")

        plot3 = surface(a, d, ddpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables",legend=false)
        savefig(plot3, "Output/Policy/Dpolicy_slice_e$(ie).png")

        plot4 = surface(a, d, dcpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Consumption",legend=false)
        savefig(plot4, "Output/Policy/Cpolicy_slice_e$(ie).png")

        plot5 = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, dpol[Int(floor(sz.nz/2)), ie, :, Int(floor(sz.nd/2))], label="e=$ie")
        end

        # Save the plot
        savefig(plot5, "Output/Policy/Dpolicy_lines.png")


        plot5a = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne#(1, 3, 7)
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, dpol[Int(floor(sz.nz/2)), ie, :, Int(floor(sz.nd/3))], label="e=$ie")
        end

        # Save the plot
        savefig(plot5a, "Output/Policy/Dpolicy_lines_dlow.png")
      
        plot5b = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne#(1, 3, 7)
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, dpol[Int(floor(sz.nz/2)), ie, :, Int(floor(2/3*sz.nd))], label="e=$ie")
        end

        # Save the plot
        savefig(plot5b, "Output/Policy/Dpolicy_lines_dhigh.png")


        # Create a new plot for the durable policy function
        plot6 = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne#(1, 3, 7)
        # Extract the policy function for the current economic state
        # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
            plot!(a, apol[Int(floor(sz.nz/2)), ie, :, Int(floor(sz.nd/2))], label="e=$ie")
        end
    

        savefig(plot6, "Output/Policy/Apolicy_lines.png")

        plot6a = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne#(1, 3, 7)
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, apol[Int(floor(sz.nz/2)), ie, :, Int(floor(sz.nd/3))], label="e=$ie")
        end

        # Save the plot
        savefig(plot6a, "Output/Policy/Apolicy_lines_dlow.png")
      
        plot6b = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne#(1, 3, 7)
            # Extract the policy function for the current economic state
            # Assuming apol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, apol[Int(floor(sz.nz/2)), ie, :, Int(floor(2/3*sz.nd))], label="e=$ie")
        end

        # Save the plot
        savefig(plot6b, "Output/Policy/Apolicy_lines_dhigh.png")


   # Choose a few e-levels to loop over (e.g., low, medium, high)

# Durable indices (used in both plots)
d_idx = [round(Int, sz.nd ÷ 2)+1, sz.nd]  # median and high durable index

for ie in 1:sz.ne
    # Asset Policy Plot
    plotA = plot(
        xlabel = "Initial Assets",
        ylabel = "Asset Policy",
        title = "Asset Policy by Durable Level",
        legend = :outerbottom,
        legend_orientation = :horizontal
    )

    for id in d_idx
        label_d = "d=$(round(d[id]; digits=2))"
        plot!(a, apol[Int(floor(sz.nz/2)), ie, :, id], label = label_d)
    end

    savefig(plotA, "Output/Policy/Apolicy_byD_fixedE$(ie).png")

    # Durable Policy Plot
    plotD = plot(
        xlabel = "Initial Assets",
        ylabel = "Durable Policy",
        title = "Durable Policy by Durable Level",
        legend = :outerbottom,
        legend_orientation = :horizontal
    )

    for id in d_idx
        label_d = "d=$(round(d[id]; digits=2))"
        plot!(a, dpol[Int(floor(sz.nz/2)), ie, :, id], label = label_d)
    end

    savefig(plotD, "Output/Policy/Dpolicy_byD_fixedE$(ie).png")
end

# Select low, medium, high indices for assets and durables
a_idx = [2, Int(floor(sz.na / 2))+1, sz.na-1]  # low, medium, high
d_idx = [2, Int(floor(sz.nd / 2))+1, sz.nd-1]  # low, medium, high

e_vals = g.ex  # actual exchange rate values

for ia in a_idx
    for id in d_idx
        # Asset Policy over exchange rates
        plotA = plot(
            e_vals,
            apol[Int(floor(sz.nz/2)), :, ia, id],
            ylims = (0, 60),
            xlabel = "Exchange Rate",
            ylabel = "Asset Policy",
            title = "a' vs e (a = $(round(a[ia], digits=2)), d = $(round(d[id], digits=2)))",
            marker = :circle,
            label = false
        )
        savefig(plotA, "Output/Policy/Apolicy_byE_a$(ia)_d$(id).png")

        # Durable Policy over exchange rates
        plotD = plot(
            e_vals,
            dpol[Int(floor(sz.nz/2)), :, ia, id],
            ylims = (0, 60),
            xlabel = "Exchange Rate",
            ylabel = "Durable Policy",
            title = "d' vs e (a = $(round(a[ia], digits=2)), d = $(round(d[id], digits=2)))",
            marker = :square,
            label = false
        )
        savefig(plotD, "Output/Policy/Dpolicy_byE_a$(ia)_d$(id).png")
    end
end



    end


end
