using Plots
default(fontfamily = "Computer Modern")

function decision_rules(answ)
    output_dir = "Output/Aggregates"
    isdir(output_dir) || mkpath(output_dir)

    # Grids
    ex = answ.g.ex
    a  = answ.g.a
    d  = answ.g.d

    # Policies & indicator
    d_pol                  = answ.pol.d
    d_adjust_pol           = answ.adjust_result.pol.d
    adjust_indicator_policy= answ.adjustment_indicator  # Bool/Bit array [ie,iy,iaa,ia,id]

    # choose a middle y and aa state, like your style
    iy0  = Int(clamp(floor(sz.ny/2), 1, sz.ny))
    iaa0 = Int(clamp(floor(sz.na/2), 1, sz.na))

    # Δd for info (adjust-only)
    d_change = d_adjust_pol .- reshape(d, (1,1,1,1,sz.nd))
    d_change_adjust = zeros(size(d_change))
    d_change_adjust[adjust_indicator_policy] .= d_change[adjust_indicator_policy]
    d_sign_adjust = sign.(d_change_adjust)

    # Fixed d: heatmap over (e × a), fixing y=iy0 and aa=iaa0
    for id in 1:sz.nd
        H = dropdims(adjust_indicator_policy[:, iy0, iaa0, :, id], dims=(2,3)) # [ie, ia]
        heatmap(a, ex, H, xlabel="Foreign assets a", ylabel="Exchange rate e",
                title="Adjust regions (fix d[$id], y=iy0, aa=mid)", color=:blues)
        savefig(joinpath(output_dir, "FixedD_Decision_Rules_d$id.png"))
    end

    # Fixed a: heatmap over (e × d), fixing y=iy0, aa=iaa0, and a=ia
    for ia in 1:sz.na
        H = dropdims(adjust_indicator_policy[:, iy0, iaa0, ia, :], dims=(2,3)) # [ie, id]
        heatmap(d, ex, H, xlabel="Durables d", ylabel="Exchange rate e",
                title="Adjust regions (fix a[$ia], y=iy0, aa=mid)", color=:blues)
        savefig(joinpath(output_dir, "FixedA_Decision_Rules_a$ia.png"))
    end

    # Size of adjustment, same slices
    cmap = cgrad([:blue, :white])
    for id in 1:sz.nd
        S = dropdims(d_change_adjust[:, iy0, iaa0, :, id], dims=(2,3))
        heatmap(a, ex, S, xlabel="Foreign assets a", ylabel="Exchange rate e",
                title="Adjustment size (fix d[$id], y=iy0, aa=mid)", color=cmap)
        savefig(joinpath(output_dir, "Size_Decision_Rules_d$id.png"))
    end
    for ia in 1:sz.na
        S = dropdims(d_change_adjust[:, iy0, iaa0, ia, :], dims=(2,3))
        heatmap(d, ex, S, xlabel="Durables d", ylabel="Exchange rate e",
                title="Adjustment size (fix a[$ia], y=iy0, aa=mid)", color=cmap)
        savefig(joinpath(output_dir, "Size_Decision_Rules_a$ia.png"))
    end

    # Sign panels
    palette = [:white, :skyblue, :blue]
    levels = [-1, 0, 1]
    for id in 1:sz.nd
        S = dropdims(d_sign_adjust[:, iy0, iaa0, :, id], dims=(2,3))
        heatmap(a, ex, S, xlabel="Foreign assets a", ylabel="Exchange rate e",
                title="Durable sign change (fix d[$id], y=iy0, aa=mid)",
                color=palette, discrete_values=levels,
                colorbar_ticks=[(l, string(l)) for l in levels])
        savefig(joinpath(output_dir, "Sign_Decision_Rules_d$id.png"))
    end
    for ia in 1:sz.na
        S = dropdims(d_sign_adjust[:, iy0, iaa0, ia, :], dims=(2,3))
        heatmap(d, ex, S, xlabel="Durables d", ylabel="Exchange rate e",
                title="Durable sign change (fix a[$ia], y=iy0, aa=mid)",
                color=palette, discrete_values=levels,
                colorbar_ticks=[(l, string(l)) for l in levels])
        savefig(joinpath(output_dir, "Sign_Decision_Rules_a$ia.png"))
    end
end
