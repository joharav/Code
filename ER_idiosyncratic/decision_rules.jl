using Plots
default(fontfamily = "Computer Modern")

function decision_rules(answ)
    output_dir = "Output/Aggregates"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Grids
    ex = answ.g.ex
    a = answ.g.a
    d = answ.g.d

    # Policies
    d_pol = answ.pol.d
    d_adjust_pol = answ.adjust_result.pol.d
    adjust_indicator_policy = answ.adjustment_indicator

    # Compute d_change and its sign
    d_change = similar(d_pol)
    for id in 1:sz.nd
        d_change[:, :, :, id] .= d_pol[:, :, :, id] .- d[id]
    end

    d_sign = sign.(d_change)

    # Keep only adjusted values
    d_sign_adjust = fill(0, size(d_change))
    d_sign_adjust[adjust_indicator_policy] .= d_sign[adjust_indicator_policy]

    d_change_adjust = fill(0.0, size(d_change))
    d_change_adjust[adjust_indicator_policy] .= d_change[adjust_indicator_policy]

    # Fixed d: loop over id, fix iz = 1
    iy = Int(floor(sz.ny / 2))  # Use the middle productivity state
    for id in 1:sz.nd
        heatmap(a, ex, adjust_indicator_policy[:, iy, :, id],
            ylabel="Exchange Rate",
            xlabel="Initial Assets",
            title="Decision Rule: Adjustment Regions, fixed d[$id]",
            color=:blues
        )
        savefig(joinpath(output_dir, "FixedD_Decision_Rules_d$id.png"))
    end

    # Fixed a: loop over ia, fix iz = 1
    for ia in 1:sz.na
        heatmap(d, ex, adjust_indicator_policy[:, iy, ia, :],
            ylabel="Exchange Rate",
            xlabel="Initial Durables",
            title="Decision Rule: Adjustment Regions, fixed a[$ia]",
            color=:blues
        )
        savefig(joinpath(output_dir, "FixedA_Decision_Rules_a$ia.png"))
    end


    cmap = cgrad([:blue, :white])

    # Plot size of adjustment (only for adjusted)
    for id in 1:sz.nd
        heatmap(a, ex, d_change_adjust[:, iy, :, id],
            ylabel="Exchange Rate",
            xlabel="Initial Assets",
            title="Decision Rule: Adjustment Size, fixed d[$id]",
            color=cmap
        )
        savefig(joinpath(output_dir, "Size_Decision_Rules_d$id.png"))
    end

    for ia in 1:sz.na
        heatmap(d, ex, d_change_adjust[:, iy, ia, :],
            ylabel="Exchange Rate",
            xlabel="Initial Durables",
            title="Decision Rule: Adjustment Size, fixed a[$ia]",
            color=cmap
        )
        savefig(joinpath(output_dir, "Size_Decision_Rules_a$ia.png"))
    end

    # Plot sign of change
    palette = [:white, :skyblue, :blue]
    levels = [-1, 0, 1]

    for id in 1:sz.nd
        heatmap(a, ex, d_sign_adjust[:, iy, :, id],
            ylabel="Exchange Rate",
            xlabel="Initial Assets",
            title="Decision Rule: Durable Sign Change, fixed d[$id]",
            color=palette,
            discrete_values=levels,
            colorbar_ticks=[(l, string(l)) for l in levels]
        )
        savefig(joinpath(output_dir, "Sign_Decision_Rules_d$id.png"))
    end

    for ia in 1:sz.na
        heatmap(d, ex, d_sign_adjust[:, iy, ia, :],
            ylabel="Exchange Rate",
            xlabel="Initial Durables",
            title="Decision Rule: Durable Sign Change, fixed a[$ia]",
            color=palette,
            discrete_values=levels,
            colorbar_ticks=[(l, string(l)) for l in levels]
        )
        savefig(joinpath(output_dir, "Sign_Decision_Rules_a$ia.png"))
    end
end
