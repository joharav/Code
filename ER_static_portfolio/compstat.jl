# ==========================================================================
# 4D MODEL: Comparative statics
# ==========================================================================

using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf
using StatsBase, Polynomials

include("durable_mod.jl")
include("collectfunctions.jl")
include("ptrue.jl")

using Main.sz, Main.settings, Main.globals, Main.dtp

commence = time()

# Moment names (must match makemoments output order)
momname = [
    "duration_mean",
    "dwealth_mean",
    "dwealth_var",
    "adj_rate",
    "dollar_share",
    "dollar_vol"
]

mom_indices = [1, 2, 3, 4, 5, 6]

# Parameter names
pname = PARAM_NAMES

# LaTeX labels
param_labels = PARAM_LABELS

mom_labels = Dict(
    "duration_mean" => "Duration (years)",
    "dwealth_mean"  => "Durable/Wealth Ratio",
    "dwealth_var"   => "Durable/Wealth Var",
    "adj_rate"      => "Adjustment Rate",
    "dollar_share"  => "Dollar Share",
    "dollar_vol"    => "Dollar Share Var"
)

# Baseline parameters
pea = ptrue(sz.nop)

nvary = 8
nparam = sz.nop
nnmom = length(momname)

# Which parameters to vary
varying_params = [2, 4, 5, 7, 9, 11, 16, 17]  # delta, sigma_e, nu, F_d, r_f, kappa, chi, F_t

# Parameter ranges
maxmin = fill(NaN, nparam, 2)
maxmin[1, :] = [0.95, 0.99]   # beta
maxmin[2, :] = [0.01, 0.10]   # delta
maxmin[3, :] = [0.70, 0.95]   # rho_e
maxmin[4, :] = [0.05, 0.30]   # sigma_e
maxmin[5, :] = [0.40, 0.70]   # nu
maxmin[6, :] = [1.0, 4.0]     # gamma
maxmin[7, :] = [0.5, 5.0]     # F_d - WIDER range
maxmin[8, :] = [0.5, 2.0]     # wage
maxmin[9, :] = [0.0, 0.05]    # r_foreign
maxmin[10,:] = [0.5, 2.0]     # p_d
maxmin[11,:] = [0.01, 0.5]    # kappa - NARROWER, LOWER range
maxmin[12,:] = [0.0, 0.4]     # tau
maxmin[13,:] = [0.5, 1.5]     # h
maxmin[14,:] = [0.7, 0.95]    # rho_y
maxmin[15,:] = [0.05, 0.30]   # sigma_y
maxmin[16,:] = [0.2, 0.8]     # chi
maxmin[17,:] = [0.1, 0.6]     # F_t

# Storage
allparams = zeros(nvary, nparam)
allmoms = zeros(nvary, nparam, nnmom)

isdir("Output/Comparative") || mkpath("Output/Comparative")

# Main loop
for iparam in varying_params
    for ivary in 1:nvary
        println("Parameter $(pname[iparam]), iteration $ivary/$nvary")
        
        # Reset to baseline
        ppp = ptrue(sz.nop)
        
        # New value along range
        lo, hi = maxmin[iparam, 1], maxmin[iparam, 2]
        glop = lo + (hi - lo) * (ivary - 1) / (nvary - 1)
        ppp[iparam] = glop
        
        # Compute moments
        try
            moms_full = momentgen(ppp)
            moms_sel = moms_full[mom_indices]
            
            allparams[ivary, iparam] = glop
            allmoms[ivary, iparam, :] .= moms_sel
        catch e
            @warn "Failed at param=$(pname[iparam]), val=$glop" exception=e
            allparams[ivary, iparam] = glop
            allmoms[ivary, iparam, :] .= NaN
        end
    end
end

degree = 2

# Plotting
for iparam in varying_params
    for imom in 1:nnmom
        pname_i = pname[iparam]
        mom_i = momname[imom]
        
        x_data = allparams[:, iparam]
        y_data = allmoms[:, iparam, imom]
        
        # Skip if all NaN
        if all(isnan.(y_data))
            continue
        end
        
        # Filter valid data
        valid = .!isnan.(y_data)
        x_valid = x_data[valid]
        y_valid = y_data[valid]
        
        if length(x_valid) < 3
            continue
        end
        
        # Raw scatter
        p_raw = scatter(x_valid, y_valid,
            legend=false,
            xlabel=pname_i,
            ylabel=mom_i,
            title="$pname_i vs $mom_i"
        )
        savefig(p_raw, "Output/Comparative/$(pname_i)_$(mom_i)_raw.png")
        
        # Polynomial fit
        try
            poly_fit = Polynomials.fit(x_valid, y_valid, min(degree, length(x_valid)-1))
            x_smooth = range(minimum(x_valid), stop=maximum(x_valid), length=100)
            y_smooth = poly_fit.(x_smooth)
            
            p_smooth = scatter(x_valid, y_valid, label="Data")
            plot!(p_smooth, x_smooth, y_smooth, 
                linewidth=2, label="Fit",
                xlabel=param_labels[pname_i],
                ylabel=mom_labels[mom_i],
                title="Effect of $(param_labels[pname_i])"
            )
            savefig(p_smooth, "Output/Comparative/$(pname_i)_$(mom_i)_smooth.png")
        catch
        end
    end
end

# Summary table
println("\n" * "="^60)
println("Comparative Statics Summary")
println("="^60)

for iparam in varying_params
    println("\n$(pname[iparam]):")
    for imom in 1:nnmom
        y = allmoms[:, iparam, imom]
        valid = .!isnan.(y)
        if any(valid)
            y_range = maximum(y[valid]) - minimum(y[valid])
            y_mean = mean(y[valid])
            println("  $(momname[imom]): range=$(round(y_range, digits=4)), mean=$(round(y_mean, digits=4))")
        end
    end
end

arret = time()
println("\nTotal time: $(round(arret - commence, digits=1)) seconds")
