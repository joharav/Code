using Plots
default(fontfamily = "Computer Modern")  # Looks like LaTeX

function plot_aggregates(simdata)
    # Extract assets and durable stock information
    assets          = simdata.a[sz.burnin-2:sz.nYears, :]  # Assuming 'a' represents asset levels, shape: (nYears, nFirms)
    durable_stock   = simdata.d[sz.burnin-2:sz.nYears, :]  # Assuming 'd' represents durable stock, shape: (nYears, nFirms)
    ex              = simdata.ex[sz.burnin-2:sz.nYears, :]

    # Sum over all firms/households for each time period
    aggregate_assets = mean(assets, dims=2)  # Sum along the columns for each time step
    aggregate_durable_stock = mean(durable_stock, dims=2)
    exchange_rate = mean(ex, dims=2)

    # Prepare time axis (assuming one data point per time period)
    time_periods = 1:size(assets, 1)
    
    # Ensure the output directory exists
    output_dir="Output/Aggregates"
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Plot aggregate assets
    p1 = plot(time_periods, aggregate_assets, label="Aggregate Assets", xlabel="Time", ylabel="Total Assets", title="Mean Assets Over Time", color=:blue, legend=false)
    savefig(p1, joinpath(output_dir, "Aggregate_Assets.png"))

    # Plot aggregate durable stock
    p2 = plot(time_periods, aggregate_durable_stock, label="Aggregate Durable Stock", xlabel="Time", ylabel="Mean Durable Stock", title="Aggregate Durable Stock Over Time", color=:green, legend=false)
    savefig(p2, joinpath(output_dir, "Aggregate_Durable_Stock.png"))

    # Plot aggregate durable stock
    p3 = plot(time_periods, exchange_rate, label="Exchange Rate", xlabel="Time", ylabel="Exchange Rate", title="Exchange Rate Over Time", color=:green, legend=false)
    savefig(p3, joinpath(output_dir, "Exchange_Rate.png"))
end

# Example call to create the plots, given `simdata` as your simulation output.
#plot_aggregates(simdata)
#plot_aggregates(simdata_shock)