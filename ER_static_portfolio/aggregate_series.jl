function plot_aggregates(simdata)
    # tail window
    T0 = max(sz.nYears-500+1, 1)
    a   = simdata.a[T0:sz.nYears, :]
    aa  = simdata.aa[T0:sz.nYears, :]
    d   = simdata.d[T0:sz.nYears, :]
    ex  = simdata.ex[T0:sz.nYears, :]

    a_eff = aa .+ ex .* a

    agg_a    = mean(a,   dims=2)
    agg_aa   = mean(aa,  dims=2)
    agg_aeff = mean(a_eff, dims=2)
    agg_d    = mean(d,   dims=2)
    ex_bar   = mean(ex,  dims=2)
    time     = 1:size(a,1)

    output_dir="Output/Aggregates"
    isdir(output_dir) || mkpath(output_dir)

    savefig(plot(time, agg_a,    xlabel="Time", ylabel="Mean a",    title="Mean foreign assets", legend=false),
            joinpath(output_dir, "Aggregate_a.png"))
    savefig(plot(time, agg_aa,   xlabel="Time", ylabel="Mean aa",   title="Mean local assets", legend=false),
            joinpath(output_dir, "Aggregate_aa.png"))
    savefig(plot(time, agg_aeff, xlabel="Time", ylabel="Mean a_eff",title="Mean effective assets", legend=false),
            joinpath(output_dir, "Aggregate_a_eff.png"))
    savefig(plot(time, agg_d,    xlabel="Time", ylabel="Mean d",    title="Mean durables", legend=false),
            joinpath(output_dir, "Aggregate_Durable_Stock.png"))
    savefig(plot(time, ex_bar,   xlabel="Time", ylabel="Exchange rate", title="Exchange rate over time", legend=false),
            joinpath(output_dir, "Exchange_Rate.png"))

    # distributions
    histogram(vec(a),  bins=50, normalize=true, xlabel="a",  ylabel="Density", title="Foreign assets", legend=false)
    savefig(joinpath(output_dir, "Assets_a_distr.png"))
    histogram(vec(aa), bins=50, normalize=true, xlabel="aa", ylabel="Density", title="Local assets", legend=false)
    savefig(joinpath(output_dir, "Assets_aa_distr.png"))
    histogram(vec(d),  bins=50, normalize=true, xlabel="d",  ylabel="Density", title="Durables", legend=false)
    savefig(joinpath(output_dir, "Durables_distr.png"))
end
