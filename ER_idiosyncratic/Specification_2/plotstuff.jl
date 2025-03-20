using Plots

function plotstuff(vee::Array{Float64, 4}, apol::Array{Int64, 4}, dpol::Array{Int64, 4},  cpol::Array{Float64, 4}, g::NamedTuple)

    a = g.a  # Assets
    d = g.d  # Durables

    # Iterate over price and exchange rate slices (fix: a, d ordering)
    for iy in 1:sz.ny 
        for ie in 1:sz.ne
            dvee_slice = vee[ie, iy, :, :]  
            dapol_slice = apol[ie, iy, :, :]
            ddpol_slice = dpol[ie, iy, :, :]
            dcpol_slice = cpol[ie, iy, :, :]

            plot1 = surface(a, d, dvee_slice, xlabel="Assets", ylabel="Durables", zlabel="Value Function")
            savefig(plot1, "Specification_2/Output/Policy/vf_slice_e$(ie)_$(iy).png")

            plot2 = surface(a, d, dapol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets")
            savefig(plot2, "Specification_2/Output/Policy/Apolicy_slice_e$(ie)_$(iy).png")

            plot3 = surface(a, d, ddpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables")
            savefig(plot3, "Specification_2/Output/Policy/Dpolicy_slice_e$(ie)_$(iy).png")

            plot4 = surface(a, d, dcpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Consumption")
            savefig(plot4, "Specification_2/Output/Policy/Cpolicy_slice_e$(ie)_$(iy).png")

        end
    end


end
