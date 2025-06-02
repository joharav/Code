using Plots
default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plotstuff(vee::Array{Float64, 3}, apol::Array{Int64, 3}, dpol::Array{Int64, 3},  cpol::Array{Float64, 3}, g::NamedTuple)

    a = g.a  # Assets
    d = g.d  # Durables

    # Iterate over price and exchange rate slices (fix: a, d ordering)
    for ie in 1:sz.ne
        dvee_slice = vee[ie, :, :]  
        dapol_slice = apol[ie, :, :]
        ddpol_slice = dpol[ie, :, :]
        dcpol_slice = cpol[ie, :, :]

        plot1 = surface(a, d, dvee_slice, xlabel="Assets", ylabel="Durables", zlabel="Value Function")
        savefig(plot1, "Specification_3/Output/Policy/vf_slice_e$(ie).png")

        plot2 = surface(a, d, dapol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets")
        savefig(plot2, "Specification_3/Output/Policy/Apolicy_slice_e$(ie).png")

        plot3 = surface(a, d, ddpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables")
        savefig(plot3, "Specification_3/Output/Policy/Dpolicy_slice_e$(ie).png")

        plot4 = surface(a, d, dcpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Consumption")
        savefig(plot4, "Specification_3/Output/Policy/Cpolicy_slice_e$(ie).png")
        plot5 = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, dpol[ie, :, Int(floor(sz.nd/2))], label="e=$ie")
        end

        # Save the plot
        savefig(plot5, "Specification_3/Output/Policy/Dpolicy_lines.png")


        plot5a = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, dpol[ie, :, Int(floor(sz.nd/3))], label="e=$ie")
        end

        # Save the plot
        savefig(plot5a, "Specification_3/Output/Policy/Dpolicy_lines_dlow.png")
      
        plot5b = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, dpol[ie, :, Int(floor(2/3*sz.nd))], label="e=$ie")
        end

        # Save the plot
        savefig(plot5b, "Specification_3/Output/Policy/Dpolicy_lines_dhigh.png")


        # Create a new plot for the durable policy function
        plot6 = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
        # Extract the policy function for the current economic state
        # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
            plot!(a, apol[ie, :, Int(floor(sz.nd/2))], label="e=$ie")
        end
    

        savefig(plot6, "Specification_3/Output/Policy/Apolicy_lines.png")

        plot6a = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
            # Extract the policy function for the current economic state
            # Assuming dpol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, apol[ie, :, Int(floor(sz.nd/3))], label="e=$ie")
        end

        # Save the plot
        savefig(plot6a, "Specification_3/Output/Policy/Apolicy_lines_dlow.png")
      
        plot6b = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy Function", legend=:topright)

        # Loop over economic states
        for ie in 1:sz.ne
            # Extract the policy function for the current economic state
            # Assuming apol[ie, :, :] gives the durable policy for all (a, d) at fixed e
                plot!(a, apol[ie, :, Int(floor(2/3*sz.nd))], label="e=$ie")
        end

        # Save the plot
        savefig(plot6b, "Specification_3/Output/Policy/Apolicy_lines_dhigh.png")

    end


end
