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
nvary = 8  
nparam = sz.nop  
nnmom = length(momname)  # Number of selected moments

# Indices for the parameters to vary: f , ft , chi (maintenance)
f_index, pd_index, chi_index = 7, 11, 9

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
total_combinations = nvary^3
allparams = zeros(total_combinations, nparam)
allmoms = zeros(total_combinations, nnmom)

# Counter to keep track of combinations
counter = 1

# Iterate over all combinations of ft, fd, and chi
for ift in 1:nvary
    for ipd in 1:nvary
        for ichi in 1:nvary
            println("Combination: f=", ift, ", ft=", ipd, ", chi=", ichi)
            
            # Reset to default parameters
            ppp = ptrue(sz.nop)
            
            # Compute new parameter values for each combination
            ppp[f_index] = (maxmin[f_index,2] - maxmin[f_index,1]) * ift / (nvary - 1.0) + maxmin[f_index,1]
            ppp[pd_index] = (maxmin[pd_index,2] - maxmin[pd_index,1]) * ipd / (nvary - 1.0) + maxmin[pd_index,1]
            ppp[chi_index] = (maxmin[chi_index,2] - maxmin[chi_index,1]) * ichi / (nvary - 1.0) + maxmin[chi_index,1]

            # Compute the moments
            moms = momentgen(ppp)

            # Store parameter values and relevant moments
            allparams[counter, :] = ppp
            allmoms[counter, :] = moms[momorder]
            counter += 1
        end
    end
end

# Plot results for each selected moment
for (i, imom) in enumerate(momorder)
    for combination in 1:counter-1
        plot_comp = Plots.plot(1:total_combinations, allmoms[:, i], 
                               xlabel="Combination", ylabel=momname[i],
                               legend=false, label=" ")

        filename = "Output/Comparative/combination_moment_$(momname[i]).png"
        savefig(plot_comp, filename)
    end
end

# Time tracking
arret = time()
println("Elapsed time in seconds = ", arret - commence)