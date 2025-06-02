using Plots
default(fontfamily = "Computer Modern")  # Looks like LaTeX

function decision_rules(answ)

    # Example of how to call the function
    output_dir = "Output/Aggregates"
    if !isdir(output_dir) 
        mkpath(output_dir)
    end

    # Define grid of state variables
    ex      = answ.g.ex
    a       = answ.g.a
    d       = answ.g.d

    # Extract policy functions
    d_pol           = answ.pol.d  # Optimal durable choice
    d_adjust_pol    = answ.adjust_result.pol.d  # Adjusted durable choice
    d_change = similar(d_pol)  # Same size as d_pol

    # Identify adjustment decision
    adjust_indicator_policy = answ.adjustment_indicator  # 1 if adjusting, 0 otherwise

    for id in 1:sz.nd
        d_change[:, :, :, id] .= d_pol[:, :, :, id] .- d[id]
    end

    d_sign = sign.(d_change)

    # Only define sign and size where there's an actual adjustment
    d_sign_adjust = fill(0, size(d_change))
    d_sign_adjust[adjust_indicator_policy] .= d_sign[adjust_indicator_policy]

    d_change_adjust = fill(0.0, size(d_change))
    d_change_adjust[adjust_indicator_policy] .= d_change[adjust_indicator_policy]




    for id in 1:sz.nd
        heatmap(a, ex, adjust_indicator_policy[:, :, :, id], 
            ylabel="Exchange Rate",
            xlabel="Initial Assets",
            title="Decision Rule: Adjustment Regions, fixed d",
            color=:blues
        )
        savefig(joinpath(output_dir, "FixedD_Decision_Rules_d$id.png"))
    end 

    for ia in 1:sz.na
        heatmap(d, ex, adjust_indicator_policy[:, :, ia, :], 
            ylabel="Exchange Rate",
            xlabel="Initial Durables",
            title="Decision Rule: Adjustment Regions, fixed a",
            color=:blues
        )
        savefig(joinpath(output_dir, "FixedA_Decision_Rules_a$ia.png"))
    end 


    mean_adjust = dropdims(mean(adjust_indicator_policy, dims=3), dims=3)
    heatmap(a, ex, mean_adjust, 
    ylabel="Exchange Rate",
    xlabel="Assets",
    color=:blues
    )
    savefig(joinpath(output_dir, "Decision_Rules_average.png"))
    cmap = cgrad([:blue, :white])

    for id in 1:sz.nd
        heatmap(a, ex, d_change_adjust[:, :, :, id], 
            ylabel="Exchange Rate",
            xlabel="Initial Assets",
            title="Decision Rule: Adjustment Regions, fixed d",
            color=cmap
        )
        savefig(joinpath(output_dir, "Size_Decision_Rules_d$id.png"))
    end 
    

    for ia in 1:sz.na
        heatmap(d, ex, d_change_adjust[:, :, ia, :], 
            ylabel="Exchange Rate",
            xlabel="Initial Durables",
            title="Decision Rule: Adjustment Regions, fixed a",
            color=cmap
        )
        savefig(joinpath(output_dir, "Size_Decision_Rules_a$ia.png"))
    end 

    palette = [:white, :skyblue]  # red = -1, white = 0, green = 1
    levels = [-1, 0]

    for id in 1:sz.nd
        heatmap(a, ex, d_sign_adjust[:, :, :, id],
            ylabel = "Exchange Rate",
            xlabel = "Initial Assets",
            title = "Decision Rule: Durable Sign Change, fixed d",
            color = palette,
            discrete_values = levels,
            colorbar_ticks = [(l, string(l)) for l in levels]
        )
        savefig(joinpath(output_dir, "Sign_Decision_Rules_d$id.png"))
    end
    
    palette = [:white, :skyblue, :blue]  # red = -1, white = 0, green = 1
    levels = [-1, 0, 1]

    for ia in 1:sz.na
        heatmap(d, ex, d_sign_adjust[:, :, ia, :],
            ylabel = "Exchange Rate",
            xlabel = "Initial Durables",
            title = "Decision Rule: Durable Sign Change, fixed a",
            color = palette,
            discrete_values = levels,
            colorbar_ticks = [(l, string(l)) for l in levels]
        )
        savefig(joinpath(output_dir, "Sign_Decision_Rules_a$ia.png"))
    end



end