function plotstuff(vee::Array{Float64, 5}, apol::Array{Int64, 5}, dpol::Array{Int64, 5}, g::NamedTuple)
    p = g.p  # Prices (11)
    w = g.w  # Wages (2)
    e = g.e  # Exchange rates (2)
    a = g.a  # Assets (51)
    d = g.d  # Durables (51)

    # Iterate over price slices (fix: a, d ordering)
    Threads.@threads for k in 1:sz.np
        dvee_slice = vee[k, :, :, :, :]  # (51, 51)
        dapol_slice = apol[k, :, :, :, :]
        ddpol_slice = dpol[k, :, :, :, :]
        #dmew_slice= mew[k, :, :, :, :]

        display(surface(a, d, dvee_slice', xlabel="Assets", ylabel="Durables", zlabel="Value Function"))
        savefig("Output/vf_slice_p$k.png")

        display(surface(a, d, dapol_slice', xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets"))
        savefig("Output/Apolicy_slice_p$k.png")

        display(surface(a, d, ddpol_slice', xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables"))
        savefig("Output/Dpolicy_slice_p$k.png")

        #display(surface(a,d,dmew_slice',color=:coolwarm,xlabel="Assets",ylabel="Durables",zlabel="μ"))
        #plotd = surface(a,d,dmew_slice',color=:coolwarm,xlabel="Assets",ylabel="Durables",zlabel="μ")
        #savefig(plotd,"Output/mew_a$k.png")
    end

    # Iterate over asset slices (fix: p, d ordering)
    Threads.@threads for k in 1:sz.npa
        vee_slice = vee[:, :, :, k, :]  # (11, 51)
        apol_slice = apol[:, :, :, k, :]
        dpol_slice = dpol[:, :, :, k, :]
       # mew_slice= mew[:, :, :, k, :]

        display(surface(p, d, vee_slice', xlabel="Prices", ylabel="Durables", zlabel="Value Function"))
        savefig("Output/vf_slice_a$k.png")

        display(surface(p, d, apol_slice', xlabel="Prices", ylabel="Durables", zlabel="Optimal Assets"))
        savefig("Output/Apolicy_slice_a$k.png")

        display(surface(p, d, dpol_slice', xlabel="Prices", ylabel="Durables", zlabel="Optimal Durables"))
        savefig("Output/Dpolicy_slice_a$k.png")

       # display(surface(p,d,mew_slice',color=:coolwarm,xlabel="Prices",ylabel="Durables",zlabel="μ"))
        #plotd = surface(p,d,mew_slice',color=:coolwarm,xlabel="Prices",ylabel="Durables",zlabel="μ")
        #savefig(plotd,"Output/mew_$k.png")
    
    end


    # Iterate over exchange rate slices (fix: p, a ordering)
    Threads.@threads for k in 1:sz.new
        vee_slice = vee[:, k, :, :, :]  # (11, 51)
        apol_slice = apol[:, k, :, :, :]
        dpol_slice = dpol[:, k, :, :, :]
        # mew_slice= mew[:, k, :, :, :]

        display(surface(p, d, vee_slice', xlabel="Prices", ylabel="Durables", zlabel="Value Function"))
        savefig("Output/vf_slice_e$k.png")

        display(surface(p, d, apol_slice', xlabel="Prices", ylabel="Durables", zlabel="Optimal Assets"))
        savefig("Output/Apolicy_slice_e$k.png")

        display(surface(p, d, dpol_slice', xlabel="Prices", ylabel="Durables", zlabel="Optimal Durables"))
        savefig("Output/Dpolicy_slice_e$k.png")

        # display(surface(p,d,mew_slice',color=:coolwarm,xlabel="Prices",ylabel="Durables",zlabel="μ"))
        #plotd = surface(p,d,mew_slice',color=:coolwarm,xlabel="Prices",ylabel="Durables",zlabel="μ")
        #savefig(plotd,"Output/mew_e$k.png")
    
    end



end
