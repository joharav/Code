using Plots
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

    # Identify adjustment decision
    adjust_indicator_policy = d_pol .!= d_adjust_pol  # 1 if adjusting, 0 otherwise

    for id in 1:sz.nd
        heatmap(a, ex, adjust_indicator_policy[:, :, id], 
            ylabel="Exchange Rate",
            xlabel="Assets",
            title="Decision Rule: Adjustment Regions",
            color=:blues
        )
        savefig(joinpath(output_dir, "Decision_Rules_d$id.png"))
    end 
    mean_adjust = dropdims(mean(adjust_indicator_policy, dims=3), dims=3)
    heatmap(a, ex, mean_adjust, 
    ylabel="Exchange Rate",
    xlabel="Assets",
    color=:blues
    )
    savefig(joinpath(output_dir, "Decision_Rules_average.png"))

   # savefig("Output/Decision_Rule_Heatmap.png")
   # scatter!(simdata.ex, simdata.a, 
   # marker=:circle, color=:red, alpha=0.3, label="Households (Simulated)"
   # )
end