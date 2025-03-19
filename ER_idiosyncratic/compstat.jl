using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2

include("durable_mod.jl")
include("collectfunctions.jl")

using Main.sz, Main.settings, Main.kst, Main.dtp 

commence = time()

# Define moments of interest
momname = ["adjustment_ratio", "mu_d_c", "mu_d_wealth", "mu_gap"]
momorder = [13, 9, 8, 10]  # Indices of the selected moments
pname = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "w", "chi", "pd", "ft", "tau", "h"]

# Get the true parameter values
pea = ptrue(sz.nop)

# Number of variations per parameter
nvary  = 8  
nparam = sz.nop  
nnmom   = length(momname)  # Number of selected moments

# Define parameters to vary
varying_params = [7, 9, 11]  

# Define parameter ranges (min, max)
maxmin = [
    0.80  0.95;  # beta (Discount factor)
    0.05  0.40;  # delta (Depreciation rate)
    0.3   0.8;   # rho_e (Persistence of exchange rate shock)
    0.2   0.80;  # sigma_e (Volatility of exchange rate shock)
    0.40  0.90;  # nu (Share parameter for nondurable consumption)
    1.00  3.00;  # gamma (Risk aversion)
    0.80  0.90;  # f (Adjustment fixed cost)
    0.50  5.00;  # w (Wage)
    0.1   0.8;   # chi (Required maintenance)
    2     8;     # pd (Price of durables)
    0.80  0.90;  # ft (Fixed cost on wage rate)
    0.10  0.60;  # tau (Tax rate)
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
        allmoms[ivary, iparam, :] = moms[momorder]  # Select relevant moments

    end
end

# Generate plots for each selected parameter and moment
for iparam in varying_params
    for (i, imom) in enumerate(momorder)
        ptitle = "Parameter: $(pname[iparam]), Moment: $(momname[i])"
        plot_comp = Plots.plot(allparams[:, iparam], allmoms[:, iparam, i], 
                               xlabel=pname[iparam], ylabel=momname[i],
                               legend=false, label=" ")

        filename = "Output/Comparative/moment_$(pname[iparam])_$(momname[i]).png"
        savefig(plot_comp, filename)
    end
end

# Time tracking
arret = time()
println("Elapsed time in seconds = ", arret - commence)
