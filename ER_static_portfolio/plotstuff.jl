using Plots
default(fontfamily = "Computer Modern")

# keep your old 4D method; add this 5D one:
function plotstuff(vee::Array{Float64,5}, apol::Array{Int64,5}, dpol::Array{Int64,5}, cpol::Array{Float64,5}, g::NamedTuple)
    a = g.a
    d = g.d
    e_vals = g.ex

    iy  = Int(floor(sz.ny / 2))
    iaa0 = Int(clamp(floor(sz.na / 2), 1, sz.na))   # fix local-asset slice

    output_dir = "Output/Policy"
    isdir(output_dir) || mkpath(output_dir)

    # 1) Surfaces over (a,d) for each e at fixed y, aa
    for ie in 1:sz.ne
        vf = permutedims(vee[ie, iy, iaa0, :, :], (2,1))   # (nd, na)
        aa = permutedims(apol[ie, iy, iaa0, :, :], (2,1))
        dd = permutedims(dpol[ie, iy, iaa0, :, :], (2,1))
        cc = permutedims(cpol[ie, iy, iaa0, :, :], (2,1))

        savefig(surface(a, d, vf, xlabel="Foreign assets a", ylabel="Durables d", zlabel="Value", legend=false),
                joinpath(output_dir, "vf_slice_e$(ie).png"))
        savefig(surface(a, d, aa, xlabel="Foreign assets a", ylabel="Durables d", zlabel="a'(a,d)", legend=false),
                joinpath(output_dir, "Apolicy_slice_e$(ie).png"))
        savefig(surface(a, d, dd, xlabel="Foreign assets a", ylabel="Durables d", zlabel="d'(a,d)", legend=false),
                joinpath(output_dir, "Dpolicy_slice_e$(ie).png"))
        savefig(surface(a, d, cc, xlabel="Foreign assets a", ylabel="Durables d", zlabel="c(a,d)", legend=false),
                joinpath(output_dir, "Cpolicy_slice_e$(ie).png"))
    end

    # 2) Durable policy by a, for a couple of d levels (fixed y, aa, e varies across figures)
    for ie in 1:sz.ne
        for d_level in [Int(floor(sz.nd / 3)), Int(floor(2 * sz.nd / 3))]
            plt = plot(xlabel="Foreign assets a", ylabel="Durable policy d'", title="Durable Policy (d index = $d_level)", legend=:topright)
            plot!(plt, a, dpol[ie, iy, iaa0, :, d_level], label="e=$ie")
            savefig(plt, joinpath(output_dir, "Dpolicy_lines_d$(d_level)_e$(ie).png"))
        end
    end

    # 3) Asset policy by a, for a couple of d levels
    for ie in 1:sz.ne
        for d_level in [Int(floor(sz.nd / 3)), Int(floor(2 * sz.nd / 3))]
            plt = plot(xlabel="Foreign assets a", ylabel="Foreign asset policy a'", title="Asset Policy (d index = $d_level)", legend=:topright)
            plot!(plt, a, apol[ie, iy, iaa0, :, d_level], label="e=$ie")
            savefig(plt, joinpath(output_dir, "Apolicy_lines_d$(d_level)_e$(ie).png"))
        end
    end

    # 4) Fix (a,d), vary e
    a_idx = [2, Int(floor(sz.na / 2))+1, sz.na-1]
    d_idx = [2, Int(floor(sz.nd / 2))+1, sz.nd-1]
    for ia in a_idx, id in d_idx
        a_val = round(a[ia], digits=2)
        d_val = round(d[id], digits=2)

        plotA = plot(e_vals, apol[:, iy, iaa0, ia, id], xlabel="Exchange rate e", ylabel="a'", title="a' vs e (a=$a_val, d=$d_val)", marker=:circle, label=false)
        plotD = plot(e_vals, dpol[:, iy, iaa0, ia, id], xlabel="Exchange rate e", ylabel="d'", title="d' vs e (a=$a_val, d=$d_val)", marker=:square, label=false)

        savefig(plotA, joinpath(output_dir, "Apolicy_byE_a$(ia)_d$(id).png"))
        savefig(plotD, joinpath(output_dir, "Dpolicy_byE_a$(ia)_d$(id).png"))
    end

    return nothing
end
