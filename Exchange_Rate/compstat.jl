using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf;

include("durable_mod.jl");
include("collectfunctions.jl");
using Main.sz, Main.settings, Main.kst, Main.dtp; 


momname = ["mu_i","v_i","mu_a","v_a","mu_c","v_c","avg_spell_length","prob_change"]; 
pname = ["beta","delta","rho","sigma","nu", "gamma","f","rr","w"];  
pea = ptrue(sz.nop);
# ============ Run stuff ===================================
nvary  = 5; # Number of variations per parameter
nparam = 9; 

# Define the subset of parameters you want to vary
varying_params = [ 2, 3, 4, 5, 6, 7]  # Example: only vary params 2, 4, 6, and 8

# Define parameter ranges (min, max) as before
maxmin = [
    0.95  0.95;  # Fixed beta  # Discount factor
    0.10  0.20;  # Fixed delta  # Depreciation rate for durable goods 1,2
    0.70  0.95;  # Fixed rho    # AR(1) persistence for durable price 1
    0.3  0.40;  # Fixed sigma  # Volatility of durable price shock
    0.40  0.80;  # Vary nu     # Share parameter for nondurable consumption
    2.00  2.00;  # Fixed gamma # Risk aversion parameter
    0.05  0.50;  # Vary f      # Adjustment fixed cost
    0.03  0.03;  # Fixed rr    # Asset return
    100   500;   # Fixed w      # Wage
]

# Initialize storage (size based on the number of varying parameters)
num_varying = length(varying_params)
allparams = zeros(nvary, nparam)  
allmoms = zeros(nvary, nparam, sz.nmom)  
used_params = zeros(nvary * num_varying, nparam)  

# Initialize counter globally
counter = 1  

# Loop only over the selected parameters
for iparam in varying_params
    for ivary in 1:nvary
        println("Varying parameter ", pname[iparam], " iteration ", ivary)
        
        # Set default parameters
        global pea = ptrue(sz.nop);
        ppp = pea;

        # Compute the new parameter value
        glop = (maxmin[iparam,2] - maxmin[iparam,1]) * ivary / (nvary - 1.0) + maxmin[iparam,1]
        ppp[iparam,1] = glop  

        # Store the used parameter values
        used_params[counter, :] = ppp'
        
        # Compute the moments
        global moms = momentgen(ppp);
        allmoms[ivary, iparam, :] = moms'; 
        allparams[ivary, iparam] = glop;

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