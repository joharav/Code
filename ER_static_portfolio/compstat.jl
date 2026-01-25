# =======================
# Comparative statics with CLEAN curves (no PARAM_* redefinition + no @L_str)
# =======================

using Random, Distributions, LinearAlgebra, Statistics, Printf
using StatsBase, Polynomials
using Plots
gr()

include("durable_mod.jl")
include("collectfunctions.jl")
include("ptrue.jl")

using Main.sz, Main.settings

commence = time()

# -----------------------
# Moments
# -----------------------
momname = [
    "duration_mean",
    "dwealth_mean",
    "dwealth_var",
    "adj_rate",
    "dollar_share",
    "dollar_vol"
]
mom_indices = 1:length(momname)
nnmom = length(momname)

mom_labels = Dict(
    "duration_mean" => "Duration",
    "dwealth_mean"  => "Durable/Wealth Ratio",
    "dwealth_var"   => "Durable/Wealth Var",
    "adj_rate"      => "Adjustment Rate",
    "dollar_share"  => "Dollar Share",
    "dollar_vol"    => "Dollar Share Var"
)

# -----------------------
# Parameters: plain strings (avoid LaTeXStrings' @L_str macro)
# -----------------------
param_labels = Dict(
    "beta"    => raw"$\beta$",
    "delta"   => raw"$\delta$",
    "rho_e"   => raw"$\rho_e$",
    "sigma_e" => raw"$\sigma_e$",
    "nu"      => raw"$\nu$",
    "gamma"   => raw"$\gamma$",
    "Fd"      => raw"$F^d$",
    "wage"    => raw"$w$",
    "r_for"   => raw"$r^{\$}$",
    "pd"      => raw"$p_d$",
    "kappa"   => raw"$\kappa$",
    "tau"     => raw"$\tau$",
    "h"       => raw"$h$",
    "rho_y"   => raw"$\rho_y$",
    "sigma_y" => raw"$\sigma_y$",
    "chi"     => raw"$\chi$",
    "Ft"      => raw"$f_t$"
)

# -----------------------
# Map NAME -> index in pea
# (edit once if your ptrue ordering differs)
# -----------------------
param_index = Dict(
    "beta"    => 1,
    "delta"   => 2,
    "rho_e"   => 3,
    "sigma_e" => 4,
    "nu"      => 5,
    "gamma"   => 6,
    "Fd"      => 7,
    "wage"    => 8,
    "r_for"   => 9,
    "pd"      => 10,
    "kappa"   => 11,
    "tau"     => 12,
    "h"       => 13,
    "rho_y"   => 14,
    "sigma_y" => 15,
    "chi"     => 16,
    "Ft"      => 17
)

# which to vary (by NAME)
varying_param_names = ["delta", "sigma_e", "nu", "Fd", "r_for", "kappa", "chi", "Ft"]

# -----------------------
# Ranges (by NAME)
# -----------------------
param_range = Dict(
    "beta"    => (0.95, 0.99),
    "delta"   => (0.01, 0.10),
    "rho_e"   => (0.70, 0.95),
    "sigma_e" => (0.05, 0.30),
    "nu"      => (0.40, 0.70),
    "gamma"   => (1.0, 4.0),
    "Fd"      => (0.5, 5.0),
    "wage"    => (0.5, 2.0),
    "r_for"   => (0.0, 0.05),
    "pd"      => (0.5, 2.0),
    "kappa"   => (0.01, 0.5),
    "tau"     => (0.0, 0.4),
    "h"       => (0.5, 1.5),
    "rho_y"   => (0.7, 0.95),
    "sigma_y" => (0.05, 0.30),
    "chi"     => (0.2, 0.8),
    "Ft"      => (0.1, 0.6)
)

# sanity checks
for nm in varying_param_names
    @assert haskey(param_labels, nm) "Missing label for '$nm' in param_labels"
    @assert haskey(param_index, nm)  "Missing index for '$nm' in param_index"
    @assert haskey(param_range, nm)  "Missing range for '$nm' in param_range"
end

# -----------------------
# Storage
# -----------------------
nvary  = 8
nparam = sz.nop

allparams = fill(NaN, nvary, nparam)
allmoms   = fill(NaN, nvary, nparam, nnmom)

outdir = "Output/Comparative"
isdir(outdir) || mkpath(outdir)

# -----------------------
# Main loop
# -----------------------
for nm in varying_param_names
    iparam = param_index[nm]
    lo, hi = param_range[nm]

    for ivary in 1:nvary
        println("Parameter $nm (idx=$iparam), iteration $ivary/$nvary")

        ppp = ptrue(sz.nop)
        glop = lo + (hi - lo) * (ivary - 1) / (nvary - 1)
        ppp[iparam] = glop

        try
            moms_full = momentgen(ppp)
            moms_sel  = moms_full[mom_indices]

            allparams[ivary, iparam] = glop
            allmoms[ivary, iparam, :] .= moms_sel
        catch e
            @warn "Failed at param=$nm (idx=$iparam), val=$glop" exception=e
            allparams[ivary, iparam] = glop
            allmoms[ivary, iparam, :] .= NaN
        end
    end
end

# -----------------------
# Smooth helper
# -----------------------
degree = 2
function poly_smooth(x::AbstractVector, y::AbstractVector; degree::Int=2, ngrid::Int=200)
    valid = isfinite.(x) .& isfinite.(y)
    x = x[valid]; y = y[valid]
    length(x) < 3 && return nothing
    deg = min(degree, length(x)-1)
    fit = Polynomials.fit(x, y, deg)
    xs  = range(minimum(x), stop=maximum(x), length=ngrid)
    ys  = fit.(xs)
    return xs, ys
end

# -----------------------
# Plotting
# -----------------------
for nm in varying_param_names
    iparam = param_index[nm]
    x_data = allparams[:, iparam]

    for imom in 1:nnmom
        mom_i  = momname[imom]
        y_data = allmoms[:, iparam, imom]

        valid = isfinite.(x_data) .& isfinite.(y_data)
        count(valid) < 3 && continue
        x_valid = x_data[valid]
        y_valid = y_data[valid]

        # RAW
        p_raw = scatter(
            x_valid, y_valid,
            legend=false,
            xlabel=param_labels[nm],
            ylabel=mom_labels[mom_i],
            title="$(param_labels[nm]) vs $(mom_labels[mom_i])"
        )
        savefig(p_raw, joinpath(outdir, "$(nm)_$(mom_i)_raw.pdf"))

        # Smooth fit
        sm = poly_smooth(x_valid, y_valid; degree=degree, ngrid=200)
        sm === nothing && continue
        xs, ys = sm

        p_smooth = scatter(
            x_valid, y_valid,
            label="Data",
            xlabel=param_labels[nm],
            ylabel=mom_labels[mom_i],
            title="Effect of $(param_labels[nm])"
        )
        plot!(p_smooth, xs, ys, linewidth=2, label="Fit")
        savefig(p_smooth, joinpath(outdir, "$(nm)_$(mom_i)_smooth.pdf"))

        # CLEAN curve only
        p_clean = plot(
            xs, ys,
            linewidth=2,
            legend=false,
            xlabel=param_labels[nm],
            ylabel=mom_labels[mom_i],
            title="Effect of $(param_labels[nm])"
        )
        savefig(p_clean, joinpath(outdir, "clean_$(nm)_$(mom_i).pdf"))
    end
end

# -----------------------
# Summary
# -----------------------
println("\n" * "="^60)
println("Comparative Statics Summary")
println("="^60)

for nm in varying_param_names
    iparam = param_index[nm]
    println("\n$nm (idx=$iparam):")
    for imom in 1:nnmom
        y = allmoms[:, iparam, imom]
        valid = isfinite.(y)
        if any(valid)
            y_range = maximum(y[valid]) - minimum(y[valid])
            y_mean  = mean(y[valid])
            println("  $(momname[imom]): range=$(round(y_range, digits=4)), mean=$(round(y_mean, digits=4))")
        end
    end
end

println("\nTotal time: $(round(time() - commence, digits=1)) seconds")
