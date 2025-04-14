using Plots
default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plot_policy_functions(policies::OrderedDict)

    output_dir = "Output/PolicyPictures/"
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # Define line styles
    linestyles = [:solid, :dash, :dot]  # Can add more styles like :dashdot, :dashdotdot if needed
    style_keys = collect(keys(policies))


    p1=plot(
        xlabel = "Current assets, a",
        ylabel = "Durable Policy, d'",
        title = "Durable Policy Function",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.d[ne_fix, :, Int(floor(sz.nd/2))+1]
            plot!(policy.g.a, avg_policy,            label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p1, joinpath(output_dir, "DurablePolicy.png"))

    p1h=plot(
        xlabel = "Current assets, a",
        ylabel = "Durable Policy, d'",
        title = "Durable Policy Function Comparison",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.d[ne_fix, :, Int(floor(sz.nd*2/3))+1]
            plot!(policy.g.a, avg_policy,            label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p1h, joinpath(output_dir, "DurablePolicy_high.png"))


    p1l=plot(
        xlabel = "Current assets, a",
        ylabel = "Durable Policy, d'",
        title = "Durable Policy Function Comparison",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.d[ne_fix, :, Int(floor(sz.nd/3))+1]
            plot!(policy.g.a, avg_policy,            label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p1l, joinpath(output_dir, "DurablePolicy_low.png"))

    p2=plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy Function",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.a[ne_fix, :, Int(floor(sz.nd/2))+1]
            plot!(policy.g.a, avg_policy,      label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p2, joinpath(output_dir, "AssetPolicy.png"))


    p2_1=plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy before adjustment",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.a[ne_fix, 1:3, Int(floor(sz.nd/2))+1]
            plot!(policy.g.a[1:3], avg_policy,      label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p2_1, joinpath(output_dir, "AssetPolicy_1.png"))

    p2_2=plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy at adjustment",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.a[ne_fix, 1:4, Int(floor(sz.nd/2))+1]
            plot!(policy.g.a[1:4], avg_policy,      label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p2_2, joinpath(output_dir, "AssetPolicy_2.png"))

    p3=plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy Function Comparison",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.a[ne_fix, :, Int(floor(sz.nd*2/3))+1]
            plot!(policy.g.a, avg_policy,      label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p3, joinpath(output_dir, "AssetPolicyComparison_high.png"))

    p4=plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy Function Comparison",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )

    for (i, label) in enumerate(style_keys)

            policy = policies[label]
            if label == "Fixed Exchange Rate"
                ne_fix = 1
            else
                ne_fix = 4
            end
            avg_policy =policy.pol.a[ne_fix, :, Int(floor(sz.nd/3))+1]
            plot!(policy.g.a, avg_policy,      label = label,
            linestyle = linestyles[mod1(i, length(linestyles))])
    end
    savefig(p4, joinpath(output_dir, "AssetPolicyComparison_low.png"))



    p5 = plot(
        xlabel = "Current assets, a",
        ylabel = "Durable Policy, d'",
        title = "Durable Policy Function by Exchange Rate State",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )
    
    nd_fix = Int(floor(sz.nd / 2)) + 1
    e_states = [1, 4, 7]
    e_labels = ["Low e", "e=1", "High e"]
    
    for (i, (ie, elabel)) in enumerate(zip(e_states, e_labels))
        label = "Baseline (Floating exchange rate)"
        policy = policies[label]
    
        durable_policy = policy.pol.d[ie, :, nd_fix]
        plot!(policy.g.a, durable_policy,
              label = elabel,
              linestyle = linestyles[mod1(i, length(linestyles))])
    end
    
    savefig(p5, joinpath(output_dir, "DurablePolicy_by_exchange_rate.png"))

    p5h = plot(
        xlabel = "Current assets, a",
        ylabel = "Durable Policy, d'",
        title = "Durable Policy Function by Exchange Rate State",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )
    
    nd_fix = Int(floor(sz.nd* 2/ 3)) + 1
    e_states = [1, 4, 7]
    e_labels = ["Low e", "e=1", "High e"]
    
    for (i, (ie, elabel)) in enumerate(zip(e_states, e_labels))
        label = "Baseline (Floating exchange rate)"
        policy = policies[label]
    
        durable_policy = policy.pol.d[ie, :, nd_fix]
        plot!(policy.g.a, durable_policy,
              label = elabel,
              linestyle = linestyles[mod1(i, length(linestyles))])
    end
    
    savefig(p5h, joinpath(output_dir, "DurablePolicy_by_exchange_rate_high.png"))


    p5l = plot(
        xlabel = "Current assets, a",
        ylabel = "Durable Policy, d'",
        title = "Durable Policy Function by Exchange Rate State",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )
    
    nd_fix = Int(floor(sz.nd* 1/ 3)) + 1
    e_states = [1, 4, 7]
    e_labels = ["Low e", "e=1", "High e"]
    
    for (i, (ie, elabel)) in enumerate(zip(e_states, e_labels))
        label = "Baseline (Floating exchange rate)"
        policy = policies[label]
    
        durable_policy = policy.pol.d[ie, :, nd_fix]
        plot!(policy.g.a, durable_policy,
              label = elabel,
              linestyle = linestyles[mod1(i, length(linestyles))])
    end
    
    savefig(p5l, joinpath(output_dir, "DurablePolicy_by_exchange_rate_low.png"))



    p6 = plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy Function by Exchange Rate State",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )
    
    nd_fix = Int(floor(sz.nd / 2)) + 1
    e_states = [1, 4, 7]
    e_labels = ["Low e", "e=1", "High e"]
    
    for (i, (ie, elabel)) in enumerate(zip(e_states, e_labels))
        label = "Baseline (Floating exchange rate)"
        policy = policies[label]
    
        asset_policy = policy.pol.a[ie, :, nd_fix]
        plot!(policy.g.a, asset_policy,
              label = elabel,
              linestyle = linestyles[mod1(i, length(linestyles))])
    end
    
    savefig(p6, joinpath(output_dir, "AssetPolicy_by_exchange_rate.png"))

    p6h = plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy Function by Exchange Rate State",
        legend = :outerbottom, legend_orientation = :horizontal,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )
    
    nd_fix = Int(floor(sz.nd* 2/ 3)) + 1
    e_states = [1, 4, 7]
    e_labels = ["Low e", "e=1", "High e"]
    
    for (i, (ie, elabel)) in enumerate(zip(e_states, e_labels))
        label = "Baseline (Floating exchange rate)"
        policy = policies[label]
    
        asset_policy = policy.pol.a[ie, :, nd_fix]
        plot!(policy.g.a, asset_policy,
              label = elabel,
              linestyle = linestyles[mod1(i, length(linestyles))])
    end
    
    savefig(p6h, joinpath(output_dir, "AssetPolicy_by_exchange_rate_high.png"))


    p6l = plot(
        xlabel = "Current assets, a",
        ylabel = "Asset Policy, a'",
        title = "Asset Policy Function by Exchange Rate State",
        legend = :outerbottom,
        xlims = (0, 50),
        ylims = (0, 60),
        size = (800, 600)
    )
    
    nd_fix = Int(floor(sz.nd* 1/ 3)) + 1
    e_states = [1, 4, 7]
    e_labels = ["Low e", "e=1", "High e"]
    
    for (i, (ie, elabel)) in enumerate(zip(e_states, e_labels))
        label = "Baseline (Floating exchange rate)"
        policy = policies[label]
    
        asset_policy = policy.pol.a[ie, :, nd_fix]
        plot!(policy.g.a, asset_policy,
              label = elabel,
              linestyle = linestyles[mod1(i, length(linestyles))])
    end
    
    savefig(p6l, joinpath(output_dir, "AssetPolicy_by_exchange_rate_low.png"))


    println("✅ Policy functions plotted")
end
