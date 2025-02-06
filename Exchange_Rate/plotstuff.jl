using Plots

function plotstuff(vee::Array{Float64, 4}, apol::Array{Int64, 4}, dpol::Array{Int64, 4}, g::NamedTuple)

    a = g.a  # Assets
    d = g.d  # Durables

    # Iterate over price and exchange rate slices (fix: a, d ordering)
    Threads.@threads for ip in 1:sz.np
        for ie in 1:sz.ne
            dvee_slice = vee[ip, ie, :, :]  # (a, d)
            dapol_slice = apol[ip, ie, :, :]
            ddpol_slice = dpol[ip, ie, :, :]

            plot1 = surface(a, d, dvee_slice, xlabel="Assets", ylabel="Durables", zlabel="Value Function")
            savefig(plot1, "Output/vf_slice_p$(ip)_e$(ie).png")

            plot2 = surface(a, d, dapol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets")
            savefig(plot2, "Output/Apolicy_slice_p$(ip)_e$(ie).png")

            plot3 = surface(a, d, ddpol_slice, xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables")
            savefig(plot3, "Output/Dpolicy_slice_p$(ip)_e$(ie).png")
        end
    end


end
