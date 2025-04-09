# Sample household indices
sample_households = 1:5  # Assume indices 1 to 5
default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plot_asset_path_between_purchases(simdata, sample_households)
    a                   = simdata.a[sz.burnin-2:sz.nYears, :]
    adjust_indicator    = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]



    plot(title="Asset Path Between Purchases", xlabel="Time", ylabel="Asset Level")
    for hh in sample_households
        path = a[:, hh]
        adjusted_path = [path[i] for i in 1:length(path) if adjust_indicator[i, hh]]
        plot!(adjusted_path, label="Household $hh")
    end
    savefig("asset_path_plot.png")
end

