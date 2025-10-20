using StatsBase
using KernelDensity
using Interpolations
using Plots



# --- helpers ---
@inline function safe_log(x; eps=1e-8)
    return log.(max.(x, eps))
end

@inline function safe_change(num, den; eps=1e-8, scale=100.0)
    den_s = max.(den, eps)
    return scale .* (num .- den_s) ./ den_s
end

# Trapezoid rule for ∫ g(x) dx from (x,f)
function trapz(x::AbstractVector{<:Real}, g::AbstractVector{<:Real})
    @assert length(x) == length(g)
    s = 0.0
    @inbounds for i in 1:length(x)-1
        dx = x[i+1] - x[i]
        s += 0.5 * (g[i+1] + g[i]) * dx
    end
    return s
end


# ---------- gap–hazard diagnostics ----------
"""
    adjustment_gaps_sim(current_d, d_star, adjust_ind)

Returns:
- gap_vec: vec(log d* - log d)
- f_x, x_values: KDE density and support of *signed* gaps
- h_x: hazard(|g|) evaluated on KDE grid via linear interpolation with flat extrapolation
- I_d_abs: ∫ |g| h(|g|) f(g) dg  (Caballero-style magnitude integral)
- mu_gap, var_gap: moments of g
- adj_rate: unconditional adjustment freq
"""
function adjustment_gaps_sim(current_d, d_star, adjustment_indicator)
    # Calculate the adjustment gaps
    gaps    = log.(d_star) .- log.(current_d)
    gap_vec = vec(gaps)

    # Flatten indicator so its length matches gap_vec
    mask = vec(adjustment_indicator .== 1)

    # Subset only adjusted cases
    adjusted_gaps = gap_vec[mask]

    # Unconditional adjustment rate
    adj_rate = mean(mask)

    # --- guard: if no adjusted obs, return safe placeholders ---
    if length(adjusted_gaps) < 2 || length(unique(adjusted_gaps)) < 2
        return adjusted_gaps, Float64[], Float64[], Float64[], NaN, NaN, NaN, adj_rate
    end

    # KDE of signed gaps (only adjusted)
    kd = kde(adjusted_gaps)
    x_values = collect(kd.x)
    f_x = kd.density

    # Hazard h(|g|) on bins of |g|
    abs_g = abs.(gap_vec)
    nb = max(20, ceil(Int, sqrt(length(abs_g))))  # sensible bin count

    gmin, gmax = minimum(abs_g), maximum(abs_g)
    if !isfinite(gmin) || !isfinite(gmax) || gmin == gmax
        return adjusted_gaps, f_x, x_values, Float64[], NaN, mean(adjusted_gaps), var(adjusted_gaps), adj_rate
    end

    edges = range(gmin, stop=gmax, length=nb+1)
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2

    # Bin index for each observation (flattened)
    bin_idx = clamp.(searchsortedlast.(Ref(edges), abs_g), 1, nb)

    # Mean adjustment rate per bin (use flattened mask)
    hazard = [mean(mask[bin_idx .== i]) for i in 1:nb]

    # Interpolate hazard(|g|) onto |x_values|
    Hin = LinearInterpolation(centers, hazard; extrapolation_bc=Flat())
    h_x = Hin.(abs.(x_values))

    # Caballero-style magnitude integral on signed support
    I_d_abs = trapz(x_values, abs.(x_values) .* h_x .* f_x)

    mu_gap  = mean(adjusted_gaps)
    var_gap = var(adjusted_gaps)

    return adjusted_gaps, f_x, x_values, h_x, I_d_abs, mu_gap, var_gap, adj_rate
end


# ---------- duration spells ----------
"""
    completed_spells(adj::AbstractVector{Bool})

"""
function completed_spells(adj::AbstractVector{Bool})
    T = length(adj)
    spells = Int[]
    run = 0
    started = false
    @inbounds for t in 1:T
        if adj[t]                 # adjustment hits → close previous run if started
            if started && run > 0
                push!(spells, run)
            end
            run = 0
            started = true       # after the first adjustment, we can count completed spells
        else
            run += 1
        end
    end
    # trailing run is right-censored → drop
    return spells
end

"""
    spells_from_panel(adj::AbstractMatrix{<:Real})

Apply `completed_spells` to each column (household) of a T×N panel.
Assumes positive values indicate adjustment.
"""
function spells_from_panel(adj::AbstractMatrix{<:Real})
    T, N = size(adj)
    out = Int[]
    @inbounds for j in 1:N
        v = @view adj[:,j]
        s = completed_spells(vec(v .> 0))
        isempty(s) || append!(out, s)
    end
    out
end