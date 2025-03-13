using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity, JLD2;

include("durable_mod.jl");
include("collectfunctions.jl");

using Main.sz, Main.settings, Main.kst, Main.dtp; 

commence = time();
momname =["adjustment_ratio"]
#momname = ["mu_d", "v_d", "mu_a", "v_a", "mu_c", "v_c", "mu_d_income", "mu_d_wealth", "mu_d_c", "mu_gap", "var_gap", "I_d", "adjustment_ratio"]
pname = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "w", "chi", "pd", "ft", "tau", "h"]
pea = ptrue(sz.nop);
# ============ Run stuff ===================================
nvary  = 8; # Number of variations per parameter
nparam = sz.nop; 

# Define the subset of parameters you want to vary
varying_params = [7,9,1, 2, 3, 4, 5, 6, 8, 10, 11, 12]  

# Define parameter ranges (min, max)
maxmin = [
    0.80  0.95;  # beta (Discount factor)
    0.05  0.40;  # delta (Depreciation rate)
    0.3  0.8;  # rho_e (Persistence of exchange rate shock)
    0.2   0.80;  # sigma_e (Volatility of exchange rate shock)
    0.40  0.90;  # nu (Share parameter for nondurable consumption)
    1.00  3.00;  # gamma (Risk aversion)
    0.10  0.80;  # f (Adjustment fixed cost)
    0.50  5.00;     # w (Wage)
    0.2   0.9;   # chi (Required maintenance)
    2   8;       # pd (Price of durables)
    0.10  0.60;  # ft (Fixed cost on wage rate)
    0.10  0.60;  # tau (Tax rate)
]

# Initialize storage (size based on the number of varying parameters)
num_varying = length(varying_params)
allparams = zeros(nvary, nparam)  
allmoms = zeros(nvary, nparam, sz.nmom)  
used_params = zeros(nvary * num_varying, nparam)  

# Initialize counter globally
counter = 1  

# Loop over parameters and moments
for iparam in varying_params
    for ivary in 1:nvary
        println("Varying parameter ", pname[iparam], " iteration ", ivary)
        
        # Set default parameters
        global pea = ptrue(sz.nop)
        ppp = pea

        # Compute the new parameter value
        glop = (maxmin[iparam,2] - maxmin[iparam,1]) * ivary / (nvary - 1.0) + maxmin[iparam,1]
        ppp[iparam] = glop  

        # Store the used parameter values
        used_params[counter, :] = ppp'

        # Compute the moments using the updated parameters
        global moms = momentgen(ppp);

        # Store the parameter values and moments
        allparams[ivary, iparam] = glop
        allmoms[ivary, iparam, :] = moms'

        global counter += 1  
    end
end

# Plot results
for iparam in varying_params
    for imom in 1:sz.nmom
        ptitle = "Parameter: $(pname[iparam]), Moment: $(momname[imom])"
        plot_comp = Plots.plot(allparams[:, iparam], allmoms[:, iparam, imom], 
                               xlabel=pname[iparam], ylabel=momname[imom],
                               legend=false, label = " ")

        filename = "Output/Comparative/moment_$(pname[iparam])_$(momname[imom]).png"
        savefig(plot_comp, filename)
    end
end
include("heatmap.jl");

arret = time();
println("elapse of time in seconds = ",arret-commence)