using Plots

default(fontfamily = "Computer Modern")  # Use LaTeX-style fonts

function plotstuff(vee::Array{Float64, 4}, apol::Array{Int64, 4}, dpol::Array{Int64, 4},  cpol::Array{Float64, 4}, g::NamedTuple)
    a = g.a  # Asset grid
    d = g.d  # Durable grid
    e_vals = g.ex  # Exchange rate grid

    iz = Int(floor(sz.nz / 2))  # Fix a productivity index

    # Ensure output directory exists
    output_dir = "Output/Policy"
    isdir(output_dir) || mkpath(output_dir)

    # -- 1. Surface plots over (a,d) for each e --
    for ie in 1:sz.ne
        # Extract 2D slices and permute to match dimensions (rows = d, cols = a)
        vf = permutedims(vee[iz, ie, :, :], (1, 2))   # (nd, na)
        aa = permutedims(apol[iz, ie, :, :], (1, 2))
        dd = permutedims(dpol[iz, ie, :, :], (1, 2))
        cc = permutedims(cpol[iz, ie, :, :], (1, 2))

        savefig(surface(a, d, vf, xlabel="Assets", ylabel="Durables", zlabel="Value Function", legend=false),
            joinpath(output_dir, "vf_slice_e$(ie).png"))
        savefig(surface(a, d, aa, xlabel="Assets", ylabel="Durables", zlabel="Optimal Assets", legend=false),
            joinpath(output_dir, "Apolicy_slice_e$(ie).png"))
        savefig(surface(a, d, dd, xlabel="Assets", ylabel="Durables", zlabel="Optimal Durables", legend=false),
            joinpath(output_dir, "Dpolicy_slice_e$(ie).png"))
        savefig(surface(a, d, cc, xlabel="Assets", ylabel="Durables", zlabel="Optimal Consumption", legend=false),
            joinpath(output_dir, "Cpolicy_slice_e$(ie).png"))
    end

     # -- 2. Durable Policy by Asset for varying d levels --
     for ie in 1:sz.ne
         for d_level in [Int(floor(sz.nd / 3)), Int(floor(2 * sz.nd / 3))]
             plt = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy (d index = $d_level)", legend=:topright)
             plot!(plt, a, dpol[iz, ie, :, d_level], label="e=$ie")
             savefig(plt, joinpath(output_dir, "Dpolicy_lines_d$(d_level)_e$(ie).png"))
         end
     end

     # -- 3. Asset Policy by Asset for varying d levels --
     for ie in 1:sz.ne
         for d_level in [Int(floor(sz.nd / 3)), Int(floor(2 * sz.nd / 3))]
             plt = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy (d index = $d_level)", legend=:topright)
             plot!(plt, a, apol[iz, ie, :, d_level], label="e=$ie")
             savefig(plt, joinpath(output_dir, "Apolicy_lines_d$(d_level)_e$(ie).png"))
         end
     end

     # -- 4. Fixed exchange rate: Vary d to show a' and d' vs a --
     d_idx = [round(Int, sz.nd ÷ 2)+1, sz.nd]  # median and high d
     for ie in 1:sz.ne
         plotA = plot(xlabel="Assets", ylabel="Asset Policy", title="Asset Policy by d, e=$(ie)", legend=:outerbottom)
         plotD = plot(xlabel="Assets", ylabel="Durable Policy", title="Durable Policy by d, e=$(ie)", legend=:outerbottom)
         for id in d_idx
             plot!(plotA, a, apol[iz, ie, :, id], label="d=$(round(d[id], digits=2))")
             plot!(plotD, a, dpol[iz, ie, :, id], label="d=$(round(d[id], digits=2))")
         end
         savefig(plotA, joinpath(output_dir, "Apolicy_byD_fixedE$(ie).png"))
         savefig(plotD, joinpath(output_dir, "Dpolicy_byD_fixedE$(ie).png"))
     end

     # -- 5. Fixed (a,d) index: Plot over exchange rate e --
     a_idx = [2, Int(floor(sz.na / 2))+1, sz.na-1]
     d_idx = [2, Int(floor(sz.nd / 2))+1, sz.nd-1]

     for ia in a_idx, id in d_idx
         a_val = round(a[ia], digits=2)
         d_val = round(d[id], digits=2)

         plotA = plot(e_vals, apol[iz, :, ia, id], xlabel="Exchange Rate", ylabel="Asset Policy", title="a' vs e (a=$a_val, d=$d_val)", marker=:circle, label=false)
         plotD = plot(e_vals, dpol[iz, :, ia, id], xlabel="Exchange Rate", ylabel="Durable Policy", title="d' vs e (a=$a_val, d=$d_val)", marker=:square, label=false)

         savefig(plotA, joinpath(output_dir, "Apolicy_byE_a$(ia)_d$(id).png"))
         savefig(plotD, joinpath(output_dir, "Dpolicy_byE_a$(ia)_d$(id).png"))
    end

    return nothing
end
