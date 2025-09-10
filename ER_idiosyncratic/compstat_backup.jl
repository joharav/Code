using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2, Polynomials, Measures

include("durable_mod.jl")
include("collectfunctions.jl")

using Main.sz,  Main.settings, Main.globals, Main.dtp; 

commence = time()

# Only the three you’re actually returning/using now:
momname = [
    "adjustment_ratio", 
    "mu_a"
]
pname = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "w", "chi", "pd", "ft", "tau", "h","rho_y","sigma_y","theta"]

default(
    fontfamily   = "Computer Modern",
    linewidth    = 3,
    framestyle   = :box, 
    grid         = :none,
    titlefontsize= 20, 
    guidefontsize= 18, 
    tickfontsize = 13,
    size         = (1200, 720),           # wider plot
    left_margin  = 14mm,                  # more room for y-label
    right_margin = 8mm,
    top_margin   = 10mm,                  # prevent title clipping
    bottom_margin= 12mm                   # room for x-label + ticks
)

param_labels = Dict(
    "beta" => "\$\\beta\$",
    "delta" => "\$\\delta\$",
    "rho_e" => "\$\\rho_e\$",
    "sigma_e" => "\$\\sigma_e\$",
    "nu" => "\$\\nu\$",
    "gamma" => "\$\\gamma\$",
    "f" => "\$F^d\$",
    "w" => "w",
    "chi" => "\$\\chi\$",
    "pd" => "\$p_d\$",
    "ft" => "\$\\phi\$",
    "tau" => "\$\\tau\$",
    "h" => "h",
    "rho_y" => "\$\\rho_y\$",
    "sigma_y" => "\$\\sigma_y\$",
    "theta" => "\$\\theta\$"

)

mom_labels = Dict(
    "adjustment_ratio" => "Adj. frequency",
    "mu_a" => "Change in average dollar assets"
)
sel_index = Dict(
    "adjustment_ratio" => 3,
    "mu_a"             => 4,
)

# Get the true parameter values
pea = ptrue(sz.nop)
moms_baseline = momentgen(pea)
base_adj = moms_baseline[sel_index["adjustment_ratio"]]
base_mua = moms_baseline[sel_index["mu_a"]]

# Number of variations per parameter
nvary  = 8 
nparam = sz.nop  
nnmom   = length(momname)  # Number of selected moments

# Define parameters to vary
varying_params = [7]   #7,4,5, 7, 9, 11, 14, 15

# Define parameter ranges (min, max)
maxmin = [
    0.80  0.95;  # beta (Discount factor)
    0.05  0.40;  # delta (Depreciation rate)
    0.3   0.9;   # rho_e (Persistence of exchange rate shock)
    0.22   0.90;  # sigma_e (Volatility of exchange rate shock)
    0.20  0.90;  # nu (Share parameter for nondurable consumption)
    0.50  3.00;  # gamma (Risk aversion)
    0.01  0.2;  # f (Adjustment fixed cost)
    0.50  5.00;  # w (Wage)
    0.1   0.9;   # chi (Required maintenance)
    2     8;     # pd (Price of durables)
    0.10  0.95;  # ft (Fixed cost on wage rate)
    0.10  0.60;  # tau (Tax rate)
    0.1   0.5    # h (Hours worked)
    0.1   0.95;  # rho_y (AR(1) persistence for idiosyncratic income)
    0.05  0.8    # sigma_y (Volatility of idiosyncratic income shock)
    0.0   1.0     # theta (New parameter for specification 2)
]





# Storage for parameter values and selected moments
allparams = zeros(nvary, nparam)  
allmoms = zeros(nvary, nparam, nnmom)  # Store only the selected moments

# Loop over selected parameters
for iparam in varying_params
    for ivary in 1:nvary
        println("Varying parameter ", pname[iparam], " iteration ", ivary)

        # Reset to default parameters
        ppp = ptrue(sz.nop)

        # Compute new parameter value
        glop = (maxmin[iparam,2] - maxmin[iparam,1]) * ivary / (nvary - 1.0) + maxmin[iparam,1]
        println("New value for ", pname[iparam], " = ", glop)
        ppp[iparam] = glop  

        # Compute the moments
        moms = momentgen(ppp)
        mvec = [moms[sel_index[name]] for name in momname]  # length 4

        # Store parameter values and selected moments
        allparams[ivary, iparam] = glop
        allmoms[ivary, iparam, :] = mvec  # Select relevant moments

    end
end
degree = 2  # Try degree 3 or 4 for better fit

# Generate plots for each selected parameter and moment
for iparam in varying_params
    for imom in 1:nnmom
        ptitle = "Parameter: $(pname[iparam]), Moment: $(momname[imom])"
        plot_comp = Plots.plot(allparams[:, iparam], allmoms[:, iparam, imom], 
                               xlabel=pname[iparam], ylabel=momname[imom],
                               legend=false, label=" ")

        filename = "Output/Comparative/check_moment_$(pname[iparam])_$(momname[imom]).pdf"
        savefig(plot_comp, filename)
        x_data = allparams[:, iparam]  # Parameter values
        y_data = allmoms[:, iparam, imom]  # Moment values

        # Fit a higher-degree polynomial
        poly_fit = Polynomials.fit(x_data, y_data, degree)

        # Generate smooth x values
        x_smooth = range(minimum(x_data), stop=maximum(x_data), length=100)
        y_smooth = poly_fit.(x_smooth)

        # Plot
        plot_comp_smooth = plot(x_data, y_data, seriestype=:scatter, label="Raw Data", xlabel=pname[iparam], ylabel=momname[imom])
        plot!(x_smooth, y_smooth, linewidth=2, label="Polynomial Fit (deg=$degree)", linestyle=:dash)

        # Save plot
        filename = "Output/Comparative/check_smoothed_moment_$(pname[iparam])_$(momname[imom]).pdf"
        savefig(plot_comp_smooth,filename)

        # Generate smoothed fit without scatter
        plot_only_smooth = plot(
        x_smooth,
        y_smooth,
        linewidth = 2,
        label = false,
        xlabel = param_labels[pname[iparam]],
        ylabel = mom_labels[momname[imom]],
        title = "Effect of $(param_labels[pname[iparam]]) on $(mom_labels[momname[imom]])"
        )

         # Add the baseline as a horizontal dashed line
         baseline_val = (momname[imom] == "adjustment_ratio") ? base_adj : base_mua
         hline!(plot_only_smooth, [baseline_val], linestyle = :dash, color = :black, label = "Baseline")
         savefig(plot_only_smooth, "Output/Comparative/check_clean_moment_$(pname[iparam])_$(momname[imom]).pdf")

    
    end
end

# Time tracking
arret = time()
println("Elapsed time in seconds = ", arret - commence)

