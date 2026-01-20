# =============================================================================
# DIAGNOSTICS / PLOTS FOR ER_static_portfolio (4D state: e,y,w,d)
# This version is consistent with your CURRENT simulation output:
#   simdata has: v,w,d,a,aa,s,ex,y,c,adjust_indicator
# and your CURRENT valfun output:
#   answ has: v, pol (merged), g, adjust_result, noadjust_result, pea
#
# Key fixes vs your pasted code:
# - NEVER call log on negative/zero: safe_log everywhere.
# - Remove references to simdata.d_adjust (you no longer store it).
# - Fix adjustment gaps: use d_star = policy d' from ADJUST regime at current state.
# - Hazard uses bins of |gap| from ALL observations, with adjustment indicator as event.
# - Plot functions use consistent directories and avoid huge findall on full panels.
# - decision_rules uses answ.pol (merged) and answ.adjust_result vs noadjust_result
#   to compute adjustment region, and uses pol_s etc for plots.
# =============================================================================

using Statistics
using StatsBase
using KernelDensity
using Interpolations
using Plots

# ---------- helpers ----------
@inline clamp_to_grid(x, g::AbstractVector{<:Real}) = min(max(x, first(g)), last(g))

@inline function safe_log(x; eps::Float64=1e-12)
    # scalar-safe
    return log(max(x, eps))
end

@inline function safe_log_vec(x::AbstractArray; eps::Float64=1e-12)
    return log.(max.(x, eps))
end

# Trapezoid rule for ∫ g(x) dx from (x,f) with x sorted
function trapz(x::AbstractVector{<:Real}, g::AbstractVector{<:Real})
    @assert length(x) == length(g)
    s = 0.0
    @inbounds for i in 1:length(x)-1
        dx = x[i+1] - x[i]
        s += 0.5 * (g[i+1] + g[i]) * dx
    end
    return s
end

# -----------------------------------------------------------------------------
# GAP–HAZARD DIAGNOSTICS (Caballero-style)
#
# Inputs should be T×N arrays (same shapes):
#   current_d : realized d_t
#   d_star    : target d*_t (from adjust policy evaluated at current state)
#   adj_ind   : adjustment_indicator (1 if adjust event at t)
#
# Returns:
#   gap_adj      : adjusted gaps (vector)
#   f_x, xgrid   : KDE density/support of adjusted gaps (signed)
#   h_x          : hazard(|g|) interpolated onto |xgrid|
#   I_d_abs      : ∫ |g| h(|g|) f(g) dg
#   mu_gap, var_gap, adj_rate
# -----------------------------------------------------------------------------
function adjustment_gaps_sim(current_d::AbstractMatrix,
                             d_star::AbstractMatrix,
                             adj_ind::AbstractMatrix)

    # Signed gaps (vectorized)
    gap_vec = vec(safe_log_vec(d_star) .- safe_log_vec(current_d))

    # Event mask
    mask = vec(adj_ind .> 0.0)
    adj_rate = mean(mask)

    # adjusted gaps only
    gap_adj = gap_vec[mask]

    # guard: no events or degenerate
    if length(gap_adj) < 5 || length(unique(gap_adj)) < 2
        return gap_adj, Float64[], Float64[], Float64[], NaN,
               (length(gap_adj) > 0 ? mean(gap_adj) : NaN),
               (length(gap_adj) > 1 ? var(gap_adj)  : NaN),
               adj_rate
    end

    # KDE of adjusted signed gaps
    kd = kde(gap_adj)
    xgrid = collect(kd.x)
    f_x   = kd.density

    # Empirical hazard as function of |gap| using bins over all observations
    abs_g = abs.(gap_vec)

    gmin, gmax = minimum(abs_g), maximum(abs_g)
    if !isfinite(gmin) || !isfinite(gmax) || gmin == gmax
        return gap_adj, f_x, xgrid, Float64[], NaN, mean(gap_adj), var(gap_adj), adj_rate
    end

    nb = max(30, ceil(Int, sqrt(length(abs_g))))
    edges = range(gmin, stop=gmax, length=nb+1)
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2

    # bin indices (1..nb)
    bin_idx = clamp.(searchsortedlast.(Ref(edges), abs_g), 1, nb)

    # hazard per bin = P(adjust | |gap| in bin)
    hazard = Vector{Float64}(undef, nb)
    @inbounds for i in 1:nb
        inbin = (bin_idx .== i)
        denom = count(inbin)
        hazard[i] = denom > 0 ? mean(mask[inbin]) : 0.0
    end

    Hin = LinearInterpolation(centers, hazard; extrapolation_bc=Flat())
    h_x = Hin.(abs.(xgrid))

    # Caballero magnitude integral (signed support)
    I_d_abs = trapz(xgrid, abs.(xgrid) .* h_x .* f_x)

    mu_gap  = mean(gap_adj)
    var_gap = var(gap_adj)

    return gap_adj, f_x, xgrid, h_x, I_d_abs, mu_gap, var_gap, adj_rate
end


# -----------------------------------------------------------------------------
# COMPLETED SPELLS (drop right-censoring)
# -----------------------------------------------------------------------------
function completed_spells(adj::AbstractVector{Bool})
    spells = Int[]
    run = 0
    started = false
    @inbounds for t in 1:length(adj)
        if adj[t]
            if started && run > 0
                push!(spells, run)
            end
            run = 0
            started = true
        else
            run += 1
        end
    end
    return spells
end

function spells_from_panel(adj::AbstractMatrix{<:Real})
    T, N = size(adj)
    out = Int[]
    @inbounds for j in 1:N
        s = completed_spells(vec(@view(adj[:,j]) .> 0))
        isempty(s) || append!(out, s)
    end
    out
end


# -----------------------------------------------------------------------------
# PLOT AGGREGATES (consistent with simdata: a, aa, d, ex)
# -----------------------------------------------------------------------------
# ---------- helpers ----------
@inline function _derive_a_aa(w::AbstractArray, s::AbstractArray, e::AbstractArray; eps=1e-12)
    # aa in pesos, a in dollars
    aa = (1 .- s) .* w
    a  = (s .* w) ./ max.(e, eps)
    return a, aa
end

# -----------------------------------------------------------------------------
# PLOT AGGREGATES (compatible with simdata: w,d,s,ex,y,c,adjust_indicator)
# -----------------------------------------------------------------------------
function plot_aggregates(simdata; tailT::Int=500)
    T = size(simdata.w, 1)
    N = size(simdata.w, 2)
    T0 = max(T - tailT + 1, 1)

    w  = simdata.w[T0:T, :]
    d  = simdata.d[T0:T, :]
    s  = simdata.s[T0:T, :]
    ex = simdata.ex[T0:T, :]

    a, aa = _derive_a_aa(w, s, ex)
    a_eff = aa .+ ex .* a   # equals w up to numerical error; still useful as a check

    agg_w    = vec(mean(w,     dims=2))
    agg_a    = vec(mean(a,     dims=2))
    agg_aa   = vec(mean(aa,    dims=2))
    agg_aeff = vec(mean(a_eff, dims=2))
    agg_d    = vec(mean(d,     dims=2))
    ex_bar   = vec(mean(ex,    dims=2))
    time     = 1:length(agg_w)

    outdir = "Output/Aggregates"
    isdir(outdir) || mkpath(outdir)

    savefig(plot(time, agg_w,    xlabel="Time", ylabel="Mean w",     title="Mean total liquid wealth (w)", legend=false),
            joinpath(outdir, "Aggregate_w.png"))
    savefig(plot(time, agg_a,    xlabel="Time", ylabel="Mean a",     title="Mean foreign assets (a)", legend=false),
            joinpath(outdir, "Aggregate_a.png"))
    savefig(plot(time, agg_aa,   xlabel="Time", ylabel="Mean aa",    title="Mean local assets (aa)", legend=false),
            joinpath(outdir, "Aggregate_aa.png"))
    savefig(plot(time, agg_aeff, xlabel="Time", ylabel="Mean a_eff", title="Mean effective assets (aa + e*a)", legend=false),
            joinpath(outdir, "Aggregate_a_eff.png"))
    savefig(plot(time, agg_d,    xlabel="Time", ylabel="Mean d",     title="Mean durables", legend=false),
            joinpath(outdir, "Aggregate_Durable_Stock.png"))
    savefig(plot(time, ex_bar,   xlabel="Time", ylabel="Exchange rate", title="Exchange rate over time", legend=false),
            joinpath(outdir, "Exchange_Rate.png"))

    savefig(histogram(vec(w),  bins=50, normalize=true, xlabel="w",  ylabel="Density", title="Total liquid wealth (w)", legend=false),
            joinpath(outdir, "Wealth_w_distr.png"))
    savefig(histogram(vec(a),  bins=50, normalize=true, xlabel="a",  ylabel="Density", title="Foreign assets (a)", legend=false),
            joinpath(outdir, "Assets_a_distr.png"))
    savefig(histogram(vec(aa), bins=50, normalize=true, xlabel="aa", ylabel="Density", title="Local assets (aa)", legend=false),
            joinpath(outdir, "Assets_aa_distr.png"))
    savefig(histogram(vec(d),  bins=50, normalize=true, xlabel="d",  ylabel="Density", title="Durables", legend=false),
            joinpath(outdir, "Durables_distr.png"))

    return nothing
end



# -----------------------------------------------------------------------------
# ADJUSTMENT TIMING + SIZE (no simdata.d_adjust in your new sim)
# Define "size" using realized change net of depreciation:
#   size_t = | d_t - (1-δ) d_{t-1} |
# and filter to periods where adjust_indicator==1.
# -----------------------------------------------------------------------------
function d_adjust_time_size(simdata, pea; minsize::Float64=1e-2)
    delta = pea[2]

    r0 = (sz.burnin - 2):sz.nYears
    r1 = (sz.burnin - 3):(sz.nYears - 1)

    d      = simdata.d[r0, :]
    d_lag  = simdata.d[r1, :]
    adj    = simdata.adjust_indicator[r0, :]

    outdir = "Output/Aggregates"
    isdir(outdir) || mkpath(outdir)

    timing = vec(sum(adj .> 0.0, dims=2))
    savefig(histogram(timing, bins=30, xlabel="Time (post-burnin)", ylabel="Count of adjustments", legend=false),
            joinpath(outdir, "Adjustment_Timing.png"))

    # net-of-depreciation "investment" measure
    d_spend = d .- (1.0 - delta) .* d_lag
    sizes = abs.(vec(d_spend[adj .> 0.0]))
    sizes = sizes[sizes .> minsize]

    savefig(histogram(sizes, bins=80, xlabel="|Δd net of depreciation|", ylabel="Frequency", legend=false),
            joinpath(outdir, "Adjustment_Size.png"))

    return nothing
end


# -----------------------------------------------------------------------------
# DISTRIBUTIONS BY (e,y) STATE
# Avoid huge findall on full matrix; use masks and vec directly.
# -----------------------------------------------------------------------------
function plot_simulated_d_and_a_by_state(simdata)
    r0 = (sz.burnin - 2):size(simdata.w,1)

    w  = simdata.w[r0, :]
    d  = simdata.d[r0, :]
    s  = simdata.s[r0, :]
    e  = simdata.ex[r0, :]
    y  = simdata.y[r0, :]

    a, aa = _derive_a_aa(w, s, e)

    outdir = "Output/Distr_State"
    isdir(outdir) || mkpath(outdir)

    unique_e = sort(unique(vec(e)))
    unique_y = sort(unique(vec(y)))

    for ei in unique_e, yi in unique_y
        m = (e .== ei) .& (y .== yi)

        d_vals  = vec(d[m])
        a_vals  = vec(a[m])
        aa_vals = vec(aa[m])

        isempty(d_vals) && continue

        savefig(histogram(d_vals, bins=50, normalize=true,
                          xlabel="Durables", ylabel="Density",
                          title="Durables | e=$(round(ei,digits=3)), y=$(round(yi,digits=3))",
                          legend=false),
                joinpath(outdir, "d_dist_e$(round(ei,digits=3))_y$(round(yi,digits=3)).png"))

        savefig(histogram(a_vals, bins=50, normalize=true,
                          xlabel="Foreign assets a", ylabel="Density",
                          title="a | e=$(round(ei,digits=3)), y=$(round(yi,digits=3))",
                          legend=false),
                joinpath(outdir, "a_dist_e$(round(ei,digits=3))_y$(round(yi,digits=3)).png"))

        savefig(histogram(aa_vals, bins=50, normalize=true,
                          xlabel="Local assets aa", ylabel="Density",
                          title="aa | e=$(round(ei,digits=3)), y=$(round(yi,digits=3))",
                          legend=false),
                joinpath(outdir, "aa_dist_e$(round(ei,digits=3))_y$(round(yi,digits=3)).png"))
    end
    return nothing
end



# -----------------------------------------------------------------------------
# DECISION RULES / POLICY VISUALIZATION (4D: e,y,w,d)
# Uses merged policy answ.pol and adjustment region from value comparison.
# -----------------------------------------------------------------------------
function decision_rules(answ)
    outdir = "Output/Policy"
    isdir(outdir) || mkpath(outdir)

    ex = answ.g.ex
    w  = answ.g.w
    d  = answ.g.d

    pol_w = answ.pol.w
    pol_d = answ.pol.d
    pol_s = answ.pol.s

    ne, ny, nw, nd = sz.ne, sz.ny, sz.nw, sz.nd

    # Adjustment region based on values (same criterion you print)
    adj = answ.adjust_result.v .>= answ.noadjust_result.v

    # Δd = d' - d
    d_change = similar(pol_d)
    @inbounds for id in 1:nd
        d_change[:, :, :, id] .= pol_d[:, :, :, id] .- d[id]
    end

    iy_mid = cld(ny, 2)
    iw_mid = cld(nw, 2)
    ie_mid = cld(ne, 2)

    # 1) Adjust region over (e × w), fixing y=mid, varying d
    for id in 1:nd
        Z = adj[:, iy_mid, :, id]  # (ne,nw)
        savefig(heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                        title="Adjust region | d=$(round(d[id],digits=3)), y=mid",
                        legend=false),
                joinpath(outdir, "AdjRegion_d$(id).png"))
    end

    # 2) Dollar share policy over (e × w)
    for id in 1:nd
        Z = pol_s[:, iy_mid, :, id]
        savefig(heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                        title="Dollar share s* | d=$(round(d[id],digits=3)), y=mid",
                        clims=(0,1), legend=true),
                joinpath(outdir, "DollarShare_d$(id).png"))
    end

    # 3) Savings policy w'
    for id in 1:nd
        Z = pol_w[:, iy_mid, :, id]
        savefig(heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                        title="Savings w' | d=$(round(d[id],digits=3)), y=mid",
                        legend=true),
                joinpath(outdir, "Savings_d$(id).png"))
    end

    # 4) Durable change Δd
    for id in 1:nd
        Z = d_change[:, iy_mid, :, id]
        maxabs = maximum(abs.(Z))
        savefig(heatmap(w, ex, Z, xlabel="Total wealth w", ylabel="Exchange rate e",
                        title="Δd = d' - d | d=$(round(d[id],digits=3)), y=mid",
                        clims=(-maxabs, maxabs), legend=true),
                joinpath(outdir, "DurableChange_d$(id).png"))
    end

    # 5) Policies vs w (fix e=mid,y=mid)
    for id in 1:nd
        p1 = plot(w, vec(pol_w[ie_mid, iy_mid, :, id]),
                  xlabel="Current wealth w", ylabel="Next wealth w'",
                  title="Savings policy | e=mid, y=mid, d=$(round(d[id],digits=3))",
                  legend=false)
        plot!(p1, w, w, linestyle=:dash, label=false)
        savefig(p1, joinpath(outdir, "SavingsVsW_d$(id).png"))

        p2 = plot(w, vec(pol_s[ie_mid, iy_mid, :, id]),
                  xlabel="Current wealth w", ylabel="Dollar share s*",
                  title="Dollar share | e=mid, y=mid, d=$(round(d[id],digits=3))",
                  ylims=(0,1), legend=false)
        savefig(p2, joinpath(outdir, "DollarShareVsW_d$(id).png"))
    end

    # 6) Dollar share vs e (fix w=mid,y=mid)
    for id in 1:nd
        p = plot(ex, vec(pol_s[:, iy_mid, iw_mid, id]),
                 xlabel="Exchange rate e", ylabel="Dollar share s*",
                 title="Dollar share vs ER | w=mid, y=mid, d=$(round(d[id],digits=3))",
                 ylims=(0,1), legend=false)
        savefig(p, joinpath(outdir, "DollarShareVsE_d$(id).png"))
    end

    # Print summary (same as your prior version)
    println("\n=== Policy Function Summary ===")
    println("Mean adjustment prob: ", round(mean(adj), digits=4))
    println("Mean dollar share: ", round(mean(pol_s), digits=4))
    println("Dollar share range: [$(round(minimum(pol_s),digits=4)), $(round(maximum(pol_s),digits=4))]")

    println("\nDollar share by exchange rate:")
    for ie in 1:ne
        mean_s = mean(pol_s[ie, :, :, :])
        println("  e=$(round(ex[ie],digits=3)): mean s = $(round(mean_s, digits=4))")
    end

    return nothing
end


