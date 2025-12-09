using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf
using StatsBase, KernelDensity, JLD2, Polynomials

include("durable_mod.jl")
include("collectfunctions.jl")

using Main.sz, Main.settings, Main.globals, Main.dtp

commence = time()

# ---------- moments: MUST match makemoments outmoms order ----------
# outmoms = [m1, m3, m4, m6, m7] =
# [duration_mean, dwealth_mean, dwealth_var, adj_rate, owner_share]
momname = [
    "duration_mean",
    "dwealth_mean",
    "dwealth_var",
    "adj_rate",
    "dollar_share",
    "dollar_vol"
]

# indices in momentgen(ppp) / makemoments (same order here)
mom_indices = [1, 2, 3, 4, 5, 6]

# ---------- parameter names: MUST match ptrue ordering ----------
pname = [
    "beta",        # 1
    "delta",       # 2
    "rho_e",       # 3
    "sigma_e",     # 4
    "nu",          # 5
    "gamma",       # 6
    "F_d",         # 7 (F^d)
    "w",           # 8
    "r_f",         # 9 (r_foreign)
    "p_d",         # 10
    "kappa_a",     # 11
    "tau",         # 12
    "h",           # 13
    "rho_y",       # 14
    "sigma_y",     # 15
    "chi",         # 16
    "F_t"          # 17
]

default(fontfamily = "Computer Modern", titlefont = font(10), guidefont = font(9))

param_labels = Dict(
    "beta"    => "\$\\beta\$",
    "delta"   => "\$\\delta\$",
    "rho_e"   => "\$\\rho_e\$",
    "sigma_e" => "\$\\sigma_e\$",
    "nu"      => "\$\\nu\$",
    "gamma"   => "\$\\gamma\$",
    "F_d"     => "\$F^d\$",
    "w"       => "w",
    "r_f"     => "r_f",
    "p_d"     => "\$p_d\$",
    "kappa_a" => "\$\\kappa_a\$",
    "tau"     => "\$\\tau\$",
    "h"       => "h",
    "rho_y"   => "\$\\rho_y\$",
    "sigma_y" => "\$\\sigma_y\$",
    "chi"     => "\$\\chi\$",
    "F_t"     => "\$F^t\$"
)

mom_labels = Dict(
    "duration_mean" => "Durable Duration",
    "dwealth_mean"  => "Durables-Wealth Ratio",
    "dwealth_var"   => "Durables-Wealth Dispersion",
    "adj_rate"      => "Adj. Ratio",
    "dollar_share"   => "Dollar Assets Share",
    "dollar_vol"    => "Dollar Assets Volatility"
)

# True parameter vector (17×1)
pea = ptrue(sz.nop)

nvary  = 8
nparam = sz.nop
nnmom  = length(momname)

# Which parameters to move along their ranges
# Using indices in the new ordering: gamma, delta, rho_e, sigma_e, F^d
varying_params = [2, 3, 4, 5, 7, 9, 11, 16,17]

# ---------- ranges ----------
# Safer: construct an empty (nop × 2) and only fill what you need
maxmin = fill(NaN, nparam, 2)

maxmin[1, :] = [0.80, 0.95]  # beta
maxmin[2, :] = [0.05, 0.40]  # delta
maxmin[3, :] = [0.30, 0.90]  # rho_e
maxmin[4, :] = [0.20, 0.60]  # sigma_e
maxmin[5, :] = [0.40, 0.90]  # nu
maxmin[6, :] = [0.50, 3.00]  # gamma
maxmin[7, :] = [0.02, 0.80]  # F^d
maxmin[8, :] = [0.50, 5.00]  # w
maxmin[9, :] = [0.0,  0.05]  # r_f (you can refine)
maxmin[10,:] = [2.00, 8.00]  # p_d
maxmin[11,:] = [0.00, 0.90]  # kappa_a
maxmin[12,:] = [0.00, 0.60]  # tau
maxmin[13,:] = [0.10, 0.50]  # h
maxmin[14,:] = [0.10, 0.95]  # rho_y
maxmin[15,:] = [0.05, 0.80]  # sigma_y
maxmin[16,:] = [0.00, 1.00]  # chi
maxmin[17,:] = [0.00, 1.00]  # F_t

# ---------- storage ----------
allparams = zeros(nvary, nparam)
allmoms   = zeros(nvary, nparam, nnmom)

isdir("Output/Comparative") || mkpath("Output/Comparative")

# ---------- main loop ----------
for iparam in varying_params
    for ivary in 1:nvary
        println("Varying parameter ", pname[iparam], " iteration ", ivary)

        # reset to baseline
        ppp = ptrue(sz.nop)

        # new value along range
        lo, hi = maxmin[iparam, 1], maxmin[iparam, 2]
        glop = (hi - lo) * ivary / (nvary - 1.0) + lo
        ppp[iparam] = glop

        # compute full moment vector (5×1)
        moms_full = momentgen(ppp)

        # select all 5 (here it's trivial, but keeps structure explicit)
        moms_sel = moms_full[mom_indices]

        allparams[ivary, iparam]      = glop
        allmoms[ivary, iparam, :] .= moms_sel
    end
end

degree = 2

# ---------- plotting ----------
for iparam in varying_params
    for imom in 1:nnmom
        pname_i = pname[iparam]
        mom_i   = momname[imom]

        x_data = allparams[:, iparam]
        y_data = allmoms[:, iparam, imom]

        # raw scatter
        p_raw = plot(
            x_data, y_data,
            seriestype = :scatter,
            legend     = false,
            xlabel     = pname_i,
            ylabel     = mom_i,
            title      = "Parameter: $pname_i, Moment: $mom_i"
        )
        savefig(p_raw, "Output/Comparative/moment_$(pname_i)_$(mom_i).png")

        # polynomial fit
        poly_fit = Polynomials.fit(x_data, y_data, degree)
        x_smooth = range(minimum(x_data), stop=maximum(x_data), length=100)
        y_smooth = poly_fit.(x_smooth)

        p_smooth = plot(
            x_data, y_data,
            seriestype = :scatter,
            label      = "Raw Data",
            xlabel     = pname_i,
            ylabel     = mom_i
        )
        plot!(
            p_smooth,
            x_smooth,
            y_smooth,
            linewidth = 2,
            label     = "Polynomial Fit (deg=$degree)",
            linestyle = :dash
        )
        savefig(p_smooth, "Output/Comparative/smoothed_moment_$(pname_i)_$(mom_i).png")

        # clean curve only
        p_clean = plot(
            x_smooth,
            y_smooth,
            linewidth = 2,
            label     = false,
            xlabel    = param_labels[pname_i],
            ylabel    = mom_labels[mom_i],
            title     = "Effect of $(param_labels[pname_i]) on $(mom_labels[mom_i])"
        )
        savefig(p_clean, "Output/Comparative/clean_moment_$(pname_i)_$(mom_i).png")
    end
end

arret = time()
println("Elapsed time in seconds = ", arret - commence)
