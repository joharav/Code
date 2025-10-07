using Plots

function decision_rules(answ)
    # output folder
    outdir = "Output/Aggregates"
    isdir(outdir) || mkpath(outdir)

    # Grids
    ex  = answ.g.ex
    a   = answ.g.a      # foreign asset grid
    aa  = answ.g.aa     # local  asset grid
    d   = answ.g.d

    ne, ny, na, nd = sz.ne, sz.ny, sz.na, sz.nd

    # Policies / flags (5-D: ie,iy,iaa,ia,id)
    d_pol  = answ.pol.d
    adj    = answ.adjustment_indicator  # BitArray{5}

    # Compute Δd = d' - d(id) at each state
    d_change = similar(d_pol)
    @inbounds for id in 1:nd
        d_change[:,:,:,:,id] .= d_pol[:,:,:,:,id] .- d[id]
    end
    d_sign = sign.(d_change)

    # Keep only adjusted states
    d_sign_adj   = zeros(Int,  size(d_sign));   d_sign_adj[adj]   .= Int.(d_sign[adj])
    d_change_adj = zeros(Float64, size(d_change)); d_change_adj[adj] .= d_change[adj]

    # Choose “middle” indices for fixed slices
    iy_mid  = cld(ny,2)
    iaa_mid = cld(na,2)
    ia_mid  = cld(na,2)
    ie_mid  = cld(ne,2)

    # -------- Fixed d: plot over (e × a), fixing y=mid, aa=mid --------
    for id in 1:nd
        Z = adj[:, iy_mid, iaa_mid, :, id]                # (ne, na)
        heatmap(a, ex, Z, xlabel="Foreign assets a", ylabel="Exchange rate e",
                title="Adjust region | d=$(round(d[id],digits=3)), aa=mid, y=mid",
                color=:blues, legend=false)
        savefig(joinpath(outdir, "AdjRegion_fixd$(id).png"))

        Zs = d_change_adj[:, iy_mid, iaa_mid, :, id]      # (ne, na)
        heatmap(a, ex, Zs, xlabel="Foreign assets a", ylabel="Exchange rate e",
                title="Adjust size Δd | d=$(round(d[id],digits=3)), aa=mid, y=mid",
                color=cgrad([:blue, :white, :red]), legend=false)
        savefig(joinpath(outdir, "AdjSize_fixd$(id).png"))

        Zsgn = d_sign_adj[:, iy_mid, iaa_mid, :, id]      # (ne, na)
        heatmap(a, ex, Zsgn, xlabel="Foreign assets a", ylabel="Exchange rate e",
                title="Sign(Δd) | d=$(round(d[id],digits=3)), aa=mid, y=mid",
                color=[:white, :skyblue, :blue], legend=false)
        savefig(joinpath(outdir, "AdjSign_fixd$(id).png"))
    end

    # -------- Fixed a: plot over (e × d), fixing y=mid, aa=mid --------
    for ia in 1:na
        Z = adj[:, iy_mid, iaa_mid, ia, :]                 # (ne, nd)
        heatmap(d, ex, Z, xlabel="Durables d", ylabel="Exchange rate e",
                title="Adjust region | a=$(round(a[ia],digits=3)), aa=mid, y=mid",
                color=:blues, legend=false)
        savefig(joinpath(outdir, "AdjRegion_fixa$(ia).png"))
    end

    # -------- Fixed (e,y,d): plane over (a × aa) --------
    for id in 1:nd
        Z = adj[ie_mid, iy_mid, :, :, id]                  # (na, na)
        # heatmap expects size(Z) == (length(ygrid), length(xgrid)); transpose to align (aa,y) vs (a,x) if desired
        heatmap(a, aa, Z', xlabel="Foreign assets a", ylabel="Local assets aa",
                title="Adjust region | e=mid, y=mid, d=$(round(d[id],digits=3))",
                color=:blues, legend=false)
        savefig(joinpath(outdir, "AdjRegion_plane_a_aa_d$(id).png"))
    end

    return nothing
end
