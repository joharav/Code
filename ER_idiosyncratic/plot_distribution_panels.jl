using StatsBase, KernelDensity, Plots, GLM, DataFrames

"""
plot_distribution_panels(simdata, pea; outdir="Output/Dist")

Creates three figures from simulated data:
1) Histogram + KDE of durable-to-wealth ratio
2) KDE of adjustment gaps (distance to target)
3) Scatter: effective wealth vs durable-to-wealth ratio (alpha-blended)

Saves:
- Output/Dist/durable_wealth_distribution.png
- Output/Dist/adjustment_gap_density.png
- Output/Dist/wealth_vs_durable_share.png
"""



function plot_distribution_panels(simdata::NamedTuple, pea::Vector{Float64}; outdir="Output/Dist")

    isdir(outdir) || mkpath(outdir)

    # Unpack constants (match your code)
    β     = pea[1]
    w     = pea[8]
    p_d   = pea[10]
    θ     = pea[16]
    r     = 1/β - 1

    # Slice post-burnin windows like your makemoments()
    a        = simdata.a[sz.burnin-3:sz.nYears-1, :]
    d        = simdata.d[sz.burnin-3:sz.nYears-1, :]
    ex       = simdata.ex[sz.burnin-3:sz.nYears-1, :]
    y        = simdata.y[sz.burnin-3:sz.nYears-1, :]
    d_next   = simdata.d[sz.burnin-2:sz.nYears, :]              # next-period for investment calc
    d_adj    = simdata.d_adjust[sz.burnin-2:sz.nYears, :]
    adj_ind  = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]

    # Effective wealth in *same currency* (your definition)
    a_eff = θ .* ex .* a .+ (1-θ) .* a

    # Durable *value* in dollars (use current stock; aligns with your ratios)
    Vd = p_d .* ex .* d

    # --- 1) Durable-to-wealth ratio distribution ---
    # wealth proxy: financial wealth carried forward (1+r)*a_eff + Vd
    # --- Durable-to-wealth ratio distribution (clean) ---
    wealth = (1 .+ r) .* a_eff .+ Vd
    ratio_d_wealth = vec(Vd) ./ vec(wealth)

    mask   = .!(isnan.(ratio_d_wealth) .| isinf.(ratio_d_wealth))
    rW     = ratio_d_wealth[mask]
    lo,hi  = quantile(rW, (0.005, 0.995))              # trim tails for nicer axes
    rW_trim = clamp.(rW, lo, hi)

    kd  = kde(rW_trim)
    med = median(rW_trim); iqr = quantile(rW_trim, 0.75)-quantile(rW_trim,0.25)

    emp_median = 0.5  # <-- replace with your EFHU number
    emp_mean   = 0.89  # example
    

    p1 = plot(kd.x, kd.density, lw=2, legend=false,
            xlabel="Durable value / total wealth (share)",
            ylabel="Density",
            title="Distribution of durable share of wealth")
    vline!([med], lw=2, ls=:dash)
    vline!([med], lw=2, ls=:dash, color=:blue, label="Model median")
    vline!([emp_median], lw=2, ls=:dashdot, color=:red, label="Data median")
    # small stats box
    txt = @sprintf "median = %.3f\nIQR = %.3f" med iqr
    annotate!((kd.x[end] - 0.1*(kd.x[end]-kd.x[1]),
            maximum(kd.density)*0.85,
            Plots.text(txt, 8, :left)))
    savefig(p1, joinpath(outdir,"durable_wealth_distribution.png"))


    # --- 2) Adjustment gap density (uses your gap routine) ---
    # Build investment 'gap' using your function (same call pattern as makemoments)
    # We only need the returned gap vector and x-grid for plotting.
    # --- 2) Adjustment gap density + hazard ---
    gap_vec, f_x, x_vals, h_x, I_d, mu_gap, var_gap, adj_ratio =
    adjustment_gaps_sim(d, d_adj, vec(adj_ind))

    gmask = isfinite.(gap_vec)
    g = gap_vec[gmask]
    kg = kde(g)

    # Bin gaps on the KDE support and compute empirical hazard per bin
    nbins   = 40
    edges   = range(minimum(kg.x), stop=maximum(kg.x), length=nbins+1)
    centers = collect((edges[1:end-1] .+ edges[2:end]) ./ 2)
    edgesv  = collect(edges)  # for searchsortedlast

    # map each observation to its bin index in 1:nbins
    binid = map(x -> clamp(searchsortedlast(edgesv, x), 1, nbins), g)

    # align adjustment indicator to the same mask
    adj_masked = vec(adj_ind)[gmask]

    haz = Array{Float64}(undef, nbins)
    for i in 1:nbins
    sel = (binid .== i)
    if any(sel)
        haz[i] = mean(adj_masked[sel])
    else
        haz[i] = NaN
    end
    end

    # simple smoothing (moving average with window=2 on each side)
    function movavg(x; k=2)
    n = length(x); y = similar(x)
    @inbounds for i in 1:n
        lo = max(1, i-k); hi = min(n, i+k)
        y[i] = mean(skipmissing(x[lo:hi]))
    end
    y
    end

    function smooth_hazard(gap_vec, adj_ind; nbins=100)
        # Build DataFrame
        df = DataFrame(gap = gap_vec, adjust = adj_ind)
    
        # Logistic regression: Pr(adjust=1 | gap)
        logit_model = glm(@formula(adjust ~ gap), df, Binomial(), LogitLink())
    
        # Prediction grid
        xgrid = range(quantile(gap_vec, 0.01), stop=quantile(gap_vec, 0.99), length=nbins)
        pred = predict(logit_model, DataFrame(gap=xgrid))
    
        return xgrid, pred
    end

    haz_s = movavg(haz; k=2)

    xgrid, pred = smooth_hazard(g, adj_masked)

    p2 = plot(kg.x, kg.density, lw=2, label="Density f(x)",
            xlabel="Gap to target (durables, log difference)",
            ylabel="Density", title="State dependence: adjustment probability")
    twinx()
    plot!(xgrid, pred, lw=2, ls=:dot, color=:red, label="Pr(adjust|gap)")

    savefig(p2, joinpath(outdir, "adjustment_gap_hazard.png"))



    # --- 3) Wealth vs durable share (thin scatter) ---
    # --- Binned relationship: durable share by wealth quantile ---
    Wv = vec(wealth)[mask]                # reuse mask used for rW (finite)
    rWv = rW

    qedges = quantile(Wv, 0:0.1:1.0)      # deciles
    qcent  = [(qedges[i]+qedges[i+1])/2 for i in 1:10]
    means  = zeros(10); p10 = similar(means); p90 = similar(means)

    for i in 1:10
        sel = (Wv .>= qedges[i]) .& (Wv .<= qedges[i+1])
        if any(sel)
            s = rWv[sel]
            means[i] = mean(s)
            p10[i]   = quantile(s, 0.10)
            p90[i]   = quantile(s, 0.90)
        else
            means[i] = NaN; p10[i]=NaN; p90[i]=NaN
        end
    end

    p3 = plot(qcent, means, lw=2, label="Mean durable share",
            xlabel="Wealth (deciles, midpoints)",
            ylabel="Durable value / wealth",
            title="Durable share rises with wealth (binned)")
    plot!(qcent, p10, lw=1, ls=:dash, label="P10")
    plot!(qcent, p90, lw=1, ls=:dash, label="P90")
    savefig(p3, joinpath(outdir, "durable_share_by_wealth_decile.png"))

end
