using Plots

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
        savefig(plot1, "Output/Policy/vf_slice_e$(ie).png")

        plot2 = surface(a, d, dapol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets")
        savefig(plot2, "Output/Policy/Apolicy_slice_e$(ie).png")

        plot3 = surface(a, d, ddpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables")
        savefig(plot3, "Output/Policy/Dpolicy_slice_e$(ie).png")

        plot4 = surface(a, d, dcpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Consumption")
        savefig(plot4, "Output/Policy/Cpolicy_slice_e$(ie).png")

    end


end
