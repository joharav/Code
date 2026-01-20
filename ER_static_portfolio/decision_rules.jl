# ==========================================================================
# 4D MODEL: Decision rules visualization
# State: (e, y, w, d)
# ==========================================================================

using Plots

function decision_rules(answ)
    outdir = "Output/Policy"
    isdir(outdir) || mkpath(outdir)
    
    # Grids
    ex = answ.g.ex
    w = answ.g.w
    d = answ.g.d
    s_grid = answ.g.s
    
    ne, ny, nw, nd = sz.ne, sz.ny, sz.nw, sz.nd
    
    # Policies (4D: ie, iy, iw, id)
    pol_w = answ.pol.w
    pol_d = answ.pol.d
    pol_s = answ.pol.s
    
    # Adjustment indicator
    adj = answ.adjust_result.v .>= answ.noadjust_result.v
    
    # Compute Δd = d' - d(id) at each state
    d_change = similar(pol_d)
    @inbounds for id in 1:nd
        d_change[:, :, :, id] .= pol_d[:, :, :, id] .- d[id]
    end
    d_sign = sign.(d_change)
    
    # Middle indices
    iy_mid = cld(ny, 2)
    iw_mid = cld(nw, 2)
    ie_mid = cld(ne, 2)
    
    # ========================================================================
    # 1. Adjustment region: plot over (e × w), fixing y=mid, d=each
    # ========================================================================
    for id in 1:nd
        Z = adj[:, iy_mid, :, id]  # (ne, nw)
        heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                title="Adjust region | d=$(round(d[id],digits=2)), y=mid",
                color=:blues, legend=false)
        savefig(joinpath(outdir, "AdjRegion_d$(id).png"))
    end
    
    # ========================================================================
    # 2. Dollar share policy: plot over (e × w), fixing y=mid, d=each
    # ========================================================================
    for id in 1:nd
        Z = pol_s[:, iy_mid, :, id]  # (ne, nw)
        heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                title="Dollar share s* | d=$(round(d[id],digits=2)), y=mid",
                color=:viridis, clims=(0, 1), legend=true)
        savefig(joinpath(outdir, "DollarShare_d$(id).png"))
    end
    
    # ========================================================================
    # 3. Savings policy w': plot over (e × w), fixing y=mid, d=each
    # ========================================================================
    for id in 1:nd
        Z = pol_w[:, iy_mid, :, id]  # (ne, nw)
        heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                title="Savings w' | d=$(round(d[id],digits=2)), y=mid",
                color=:viridis, legend=true)
        savefig(joinpath(outdir, "Savings_d$(id).png"))
    end
    
    # ========================================================================
    # 4. Durable adjustment size: plot over (e × w)
    # ========================================================================
    for id in 1:nd
        Z = d_change[:, iy_mid, :, id]  # (ne, nw)
        maxabs = maximum(abs.(Z))
        heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                title="Δd = d' - d | d=$(round(d[id],digits=2)), y=mid",
                color=cgrad([:blue, :white, :red]), 
                clims=(-maxabs, maxabs), legend=true)
        savefig(joinpath(outdir, "DurableChange_d$(id).png"))
    end
    
    # ========================================================================
    # 5. Policy functions vs wealth (fixing e=mid, y=mid)
    # ========================================================================
    for id in 1:nd
        # Savings
        plot(w, pol_w[ie_mid, iy_mid, :, id], label="w'",
             xlabel="Current wealth w", ylabel="Next period wealth w'",
             title="Savings policy | e=mid, y=mid, d=$(round(d[id],digits=2))")
        plot!(w, w, linestyle=:dash, label="45°", color=:gray)
        savefig(joinpath(outdir, "SavingsVsW_d$(id).png"))
        
        # Dollar share
        plot(w, pol_s[ie_mid, iy_mid, :, id], label="s*",
             xlabel="Current wealth w", ylabel="Dollar share",
             title="Dollar share | e=mid, y=mid, d=$(round(d[id],digits=2))",
             ylims=(0, 1))
        savefig(joinpath(outdir, "DollarShareVsW_d$(id).png"))
    end
    
    # ========================================================================
    # 6. Dollar share vs exchange rate (fixing w=mid, y=mid)
    # ========================================================================
    for id in 1:nd
        plot(ex, pol_s[:, iy_mid, iw_mid, id], label="s*",
             xlabel="Exchange rate e", ylabel="Dollar share",
             title="Dollar share vs ER | w=mid, y=mid, d=$(round(d[id],digits=2))",
             ylims=(0, 1))
        savefig(joinpath(outdir, "DollarShareVsE_d$(id).png"))
    end
    
    # ========================================================================
    # 7. Summary statistics
    # ========================================================================
    println("\n=== Policy Function Summary ===")
    println("Mean adjustment prob: ", round(mean(adj), digits=4))
    println("Mean dollar share: ", round(mean(pol_s), digits=4))
    println("Dollar share range: [$(round(minimum(pol_s),digits=4)), $(round(maximum(pol_s),digits=4))]")
    
    # Dollar share by ER
    println("\nDollar share by exchange rate:")
    for ie in 1:ne
        mean_s = mean(pol_s[ie, :, :, :])
        println("  e=$(round(ex[ie],digits=2)): mean s = $(round(mean_s, digits=4))")
    end
    
    return nothing
end
