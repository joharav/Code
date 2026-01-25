# =======================
# 1) EPS saving (Dist plots)
# =======================
using StatsBase, KernelDensity, Plots

# IMPORTANT for EPS:
# - Use a backend that supports vector output (GR usually works; PGFPlotsX is also great for LaTeX).
# - If EPS fails on your system, switch backend: `using PGFPlotsX; pgfplotsx()`
gr()

_fin(v) = v[isfinite.(v)]

function plot_distribution_panels(simdata::NamedTuple, pea::Vector{Float64}; outdir="Output/Dist")
    isdir(outdir) || mkpath(outdir)

    pd = pea[10]

    # Post-burn-in window (match your makemoments)
    r0 = (sz.burnin - 2):sz.nYears

    d  = simdata.d[r0, :]
    ex = simdata.ex[r0, :]

    # liquid wealth in pesos (preferred): aa + e*a
    if hasproperty(simdata, :aa) && hasproperty(simdata, :a)
        aa = simdata.aa[r0, :]
        a  = simdata.a[r0, :]
        w_liq = aa .+ ex .* a
    else
        w_liq = simdata.w[r0, :]  # fallback
    end

    # durable value in pesos
    Vd = pd .* ex .* d

    # total wealth in pesos
    Wtot = Vd .+ w_liq

    ratio = vec(Vd) ./ max.(vec(Wtot), 1e-12)
    ratio = _fin(ratio)

    # trim tails for stable KDE axes
    lo, hi = quantile(ratio, (0.005, 0.995))
    ratio_t = clamp.(ratio, lo, hi)

    kd = kde(ratio_t)

    p1 = plot(kd.x, kd.density, lw=2, legend=false,
              xlabel="Durable value / (durable value + liquid wealth)",
              ylabel="Density",
              title="Durable share of wealth (simulation)")

    # Save EPS explicitly
    savefig(p1, joinpath(outdir, "durable_wealth_distribution.pdf"))

    # ---- binned relationship: durable share by wealth decile ----
    Wv = vec(Wtot)
    rv = vec(Vd) ./ max.(Wv, 1e-12)

    m = isfinite.(Wv) .& isfinite.(rv)
    Wv = Wv[m]; rv = rv[m]

    qedges = quantile(Wv, 0:0.1:1.0)
    qcent  = [(qedges[i] + qedges[i+1]) / 2 for i in 1:10]
    means  = fill(NaN, 10)
    p10    = fill(NaN, 10)
    p90    = fill(NaN, 10)

    for i in 1:10
        sel = (Wv .>= qedges[i]) .& (Wv .< qedges[i+1])  # strict upper bound avoids double-count
        if any(sel)
            s = rv[sel]
            means[i] = mean(s)
            p10[i]   = quantile(s, 0.10)
            p90[i]   = quantile(s, 0.90)
        end
    end

    p2 = plot(qcent, means, lw=2, label="Mean",
              xlabel="Total wealth (decile midpoints, pesos)",
              ylabel="Durable share",
              title="Durable share by wealth decile")
    plot!(qcent, p10, lw=1, ls=:dash, label="P10")
    plot!(qcent, p90, lw=1, ls=:dash, label="P90")

    savefig(p2, joinpath(outdir, "durable_share_by_wealth_decile.pdf"))

    return nothing
end
