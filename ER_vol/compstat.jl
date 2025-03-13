using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2

include("durable_mod.jl")
include("collectfunctions.jl")

using Main.sz, Main.settings, Main.kst, Main.dtp 

commence = time()

# Define moment of interest
momname = ["adjustment_ratio"]
pname = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "w", "chi", "pd", "ft", "tau", "h"]

# Get the true parameter values
pea = ptrue(sz.nop)

# Number of variations per parameter
nvary  = 8  
nparam = sz.nop  

# Define parameters to vary
varying_params = [7, 9, 1, 2, 3, 4, 5, 6, 8, 10, 11, 12]  

# Define parameter ranges (min, max)
maxmin = [
    0.80  0.95;  # beta (Discount factor)
    0.05  0.40;  # delta (Depreciation rate)
    0.3  0.8;    # rho_e (Persistence of exchange rate shock)
    0.2  0.80;   # sigma_e (Volatility of exchange rate shock)
    0.40  0.90;  # nu (Share parameter for nondurable consumption)
    1.00  3.00;  # gamma (Risk aversion)
    0.10  0.80;  # f (Adjustment fixed cost)
    0.50  5.00;  # w (Wage)
    0.2   0.9;   # chi (Required maintenance)
    2     8;     # pd (Price of durables)
    0.10  0.60;  # ft (Fixed cost on wage rate)
    0.10  0.60;  # tau (Tax rate)
]

# Storage for parameter values and the adjustment ratio
allparams = zeros(nvary, nparam)  
adjustment_ratios = zeros(nvary, nparam)  # Only storing this moment

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

        # Store only the adjustment ratio
        allparams[ivary, iparam] = glop
        adjustment_ratios[ivary, iparam] = moms[1]  # Assuming adjustment_ratio is the first moment

    end
end

# Plot results for the adjustment ratio
for iparam in varying_params
    ptitle = "Parameter: $(pname[iparam]), Moment: Adjustment Ratio"
    plot_comp = Plots.plot(allparams[:, iparam], adjustment_ratios[:, iparam], 
                           xlabel=pname[iparam], ylabel="Adjustment Ratio",
                           legend=false, label=" ")

    filename = "Output/Comparative/moment_$(pname[iparam])_adjustment_ratio.png"
    savefig(plot_comp, filename)
end

# Time tracking
arret = time()
println("Elapsed time in seconds = ", arret - commence)
