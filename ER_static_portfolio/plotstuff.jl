using Plots
default(fontfamily = "Computer Modern")

function plotstuff4D(
    vee::Array{Float64,4},
    pol_w::Array{Float64,4},
    pol_d::Array{Float64,4},
    pol_s::Array{Float64,4},
    pol_c::Array{Float64,4},
    g::NamedTuple;
    outdir::String = "Output/Policy"
)
    w = g.w
    d = g.d
    e_vals = g.ex

    iy = cld(sz.ny, 2)  # middle income state

    isdir(outdir) || mkpath(outdir)

    # ------------------------------------------------------------
    # 1) Surfaces over (w,d) for each e at fixed y
    # ------------------------------------------------------------
    for ie in 1:sz.ne
        vf = permutedims(vee[ie, iy, :, :], (2,1))     # (nd, nw) for surface(w,d, Z)
        ww = permutedims(pol_w[ie, iy, :, :], (2,1))
        dd = permutedims(pol_d[ie, iy, :, :], (2,1))
        ss = permutedims(pol_s[ie, iy, :, :], (2,1))
        cc = permutedims(pol_c[ie, iy, :, :], (2,1))

        savefig(surface(w, d, vf, xlabel="Total wealth w", ylabel="Durables d", zlabel="Value",
                        title="V(w,d) | e index=$ie, y=mid", legend=false),
                joinpath(outdir, "vf_slice_e$(ie).png"))

        savefig(surface(w, d, ww, xlabel="Total wealth w", ylabel="Durables d", zlabel="w'(w,d)",
                        title="Savings policy w' | e index=$ie, y=mid", legend=false),
                joinpath(outdir, "Wpolicy_slice_e$(ie).png"))

        savefig(surface(w, d, dd, xlabel="Total wealth w", ylabel="Durables d", zlabel="d'(w,d)",
                        title="Durable policy d' | e index=$ie, y=mid", legend=false),
                joinpath(outdir, "Dpolicy_slice_e$(ie).png"))

        savefig(surface(w, d, ss, xlabel="Total wealth w", ylabel="Durables d", zlabel="s(w,d)",
                        title="Dollar share s | e index=$ie, y=mid", legend=false),
                joinpath(outdir, "Spolicy_slice_e$(ie).png"))

        savefig(surface(w, d, cc, xlabel="Total wealth w", ylabel="Durables d", zlabel="c(w,d)",
                        title="Consumption c | e index=$ie, y=mid", legend=false),
                joinpath(outdir, "Cpolicy_slice_e$(ie).png"))
    end

    # ------------------------------------------------------------
    # 2) Durable policy and dollar share vs w at a couple d levels
    # ------------------------------------------------------------
    d_levels = unique(clamp.([floor(Int, sz.nd/3), floor(Int, 2*sz.nd/3)], 1, sz.nd))

    for ie in 1:sz.ne
        for id in d_levels
            d_val = round(d[id], digits=3)

            pD = plot(xlabel="Total wealth w", ylabel="d'",
                      title="Durable policy d'(w) | e index=$ie, y=mid, d=$(d_val)",
                      legend=false)
            plot!(pD, w, pol_d[ie, iy, :, id])
            savefig(pD, joinpath(outdir, "Dpolicy_lines_d$(id)_e$(ie).png"))

            pS = plot(xlabel="Total wealth w", ylabel="s",
                      title="Dollar share s(w) | e index=$ie, y=mid, d=$(d_val)",
                      ylims=(0,1), legend=false)
            plot!(pS, w, pol_s[ie, iy, :, id])
            savefig(pS, joinpath(outdir, "Spolicy_lines_d$(id)_e$(ie).png"))

            pW = plot(xlabel="Total wealth w", ylabel="w'",
                      title="Savings policy w'(w) | e index=$ie, y=mid, d=$(d_val)",
                      legend=false)
            plot!(pW, w, pol_w[ie, iy, :, id])
            plot!(pW, w, w, ls=:dash)  # 45-degree
            savefig(pW, joinpath(outdir, "Wpolicy_lines_d$(id)_e$(ie).png"))
        end
    end

    # ------------------------------------------------------------
    # 3) Fix (w,d), vary e: policies vs exchange rate
    # ------------------------------------------------------------
    w_idx = unique(clamp.([2, floor(Int, sz.nw/2)+1, sz.nw-1], 1, sz.nw))
    d_idx = unique(clamp.([2, floor(Int, sz.nd/2)+1, sz.nd-1], 1, sz.nd))

    for iw in w_idx, id in d_idx
        w_val = round(w[iw], digits=3)
        d_val = round(d[id], digits=3)

        pW = plot(e_vals, pol_w[:, iy, iw, id],
                  xlabel="Exchange rate e", ylabel="w'",
                  title="w' vs e | y=mid, w=$(w_val), d=$(d_val)",
                  marker=:circle, legend=false)
        savefig(pW, joinpath(outdir, "Wpolicy_byE_w$(iw)_d$(id).png"))

        pD = plot(e_vals, pol_d[:, iy, iw, id],
                  xlabel="Exchange rate e", ylabel="d'",
                  title="d' vs e | y=mid, w=$(w_val), d=$(d_val)",
                  marker=:square, legend=false)
        savefig(pD, joinpath(outdir, "Dpolicy_byE_w$(iw)_d$(id).png"))

        pS = plot(e_vals, pol_s[:, iy, iw, id],
                  xlabel="Exchange rate e", ylabel="s",
                  title="s vs e | y=mid, w=$(w_val), d=$(d_val)",
                  marker=:diamond, ylims=(0,1), legend=false)
        savefig(pS, joinpath(outdir, "Spolicy_byE_w$(iw)_d$(id).png"))

        pC = plot(e_vals, pol_c[:, iy, iw, id],
                  xlabel="Exchange rate e", ylabel="c",
                  title="c vs e | y=mid, w=$(w_val), d=$(d_val)",
                  marker=:utriangle, legend=false)
        savefig(pC, joinpath(outdir, "Cpolicy_byE_w$(iw)_d$(id).png"))
    end

    return nothing
end
