using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2, Polynomials

include("durable_mod.jl")
include("collectfunctions.jl")

using Main.sz,  Main.settings, Main.globals, Main.dtp; 

commence = time()

# # Define moments of interest
# momname = [
#     "mu_d",           # 1
#     "var_d",          # 2
#     "mu_a",           # 3
#     "var_a",          # 4
#     "mu_c",           # 5
#     "var_c",          # 6
#     "mu_d_wealth",    # 7
#     "mu_d_c",         # 8
#     "mu_gap",         # 9
#     "var_gap",        # 10
#     "cev",            # 11  <-- new addition
#     "adjustment_ratio", # 12
#     "IQR_d_wealth",   # 13
#     "IQR_d_c",        # 14
#     "p90_10_d_wealth",# 15
#     "p90_10_d_c"      # 16
# ]

# Only the three you’re actually returning/using now:
momname = [
    "d_dispersion",
    "mu_d_wealth",
    "adjustment_ratio"
]
pname = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "w", "chi", "pd", "ft", "tau", "h","rho_y","sigma_y","theta"]

default(fontfamily = "Computer Modern", titlefont = font(10), guidefont = font(9))  # Already done!

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

# Map raw moment names to prettier labels
# mom_labels = Dict(
#     "mu_d" => "Mean Durable rate",
#     "var_d" => "Var Durable",
#     "mu_a" => "Mean Assets rate",
#     "var_a" => "Var Assets",
#     "mu_c" => "Mean Consumption",
#     "var_c" => "Var Consumption",
#     "mu_d_wealth" => "Durables-Wealth Ratio",
#     "mu_d_c" => "Durables-Cons. Ratio",
#     "mu_gap" => "Mean Gap",
#     "var_gap" => "Var Gap",
#     "cev" => "Welfare change",
#     "adjustment_ratio" => "Adj. Ratio",
#     "IQR_d_wealth" => "Interquartile ratio durable-wealth",
#     "IQR_d_c" => "Interquartile ratio durable-consumption",
#     "p90_10_d_wealth" => "90-10 ratio durable-wealth",
#     "p90_10_d_c" => "90-10 ratio durable-consumption"

# )
mom_labels = Dict(
   "d_dispersion" => "Durable Dispersion",
    "adjustment_ratio" => "Adj. Ratio",
    "mu_d_wealth" => "Durables-Wealth Ratio"
)

# Get the true parameter values
pea = ptrue(sz.nop)

# Number of variations per parameter
nvary  = 8 
nparam = sz.nop  
nnmom   = length(momname)  # Number of selected moments

# Define parameters to vary
varying_params = [6]   #7,4,5, 7, 9, 11, 14, 15

# Define parameter ranges (min, max)
maxmin = [
    0.80  0.95;  # beta (Discount factor)
    0.05  0.40;  # delta (Depreciation rate)
    0.3   0.9;   # rho_e (Persistence of exchange rate shock)
    0.2   0.60;  # sigma_e (Volatility of exchange rate shock)
    0.20  0.90;  # nu (Share parameter for nondurable consumption)
    0.50  3.00;  # gamma (Risk aversion)
    0.05  0.8;  # f (Adjustment fixed cost)
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
        ppp[iparam] = glop  

        # Compute the moments
        moms = momentgen(ppp)

        # Store parameter values and selected moments
        allparams[ivary, iparam] = glop
        allmoms[ivary, iparam, :] = moms'  # Select relevant moments

    end
end
degree = 3  # Try degree 3 or 4 for better fit

# Generate plots for each selected parameter and moment
for iparam in varying_params
    for imom in 1:nnmom
        ptitle = "Parameter: $(pname[iparam]), Moment: $(momname[imom])"
        plot_comp = Plots.plot(allparams[:, iparam], allmoms[:, iparam, imom], 
                               xlabel=pname[iparam], ylabel=momname[imom],
                               legend=false, label=" ")

        filename = "Output/Comparative/moment_$(pname[iparam])_$(momname[imom]).png"
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
        filename = "Output/Comparative/smoothed_moment_$(pname[iparam])_$(momname[imom]).png"
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

        savefig(plot_only_smooth, "Output/Comparative/clean_moment_$(pname[iparam])_$(momname[imom]).png")

    end
end

# Time tracking
arret = time()
println("Elapsed time in seconds = ", arret - commence)

