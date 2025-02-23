using Random, Distributions, LinearAlgebra, Plots, Statistics, Printf, StatsBase, KernelDensity;

include("durable_mod.jl");
include("collectfunctions.jl");
include("plotgaps_comp.jl");

using Main.sz, Main.settings, Main.kst, Main.dtp; 

commence = time();

momname = ["mu_d", "v_d", "mu_a", "v_a", "mu_c", "v_c", "ratio_d_income", "ratio_d_wealth", "ratio_d_consumption", "mu_gap", "var_gap", "mu_hx", "var_hx", "I_d"]

pname = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "w", "chi", "pd"]
pea = ptrue(sz.nop);
# ============ Run stuff ===================================
nvary  = 5; # Number of variations per parameter
nparam = sz.nop; 

# Define the subset of parameters you want to vary
varying_params = [7, 9, 3, 4]  # Example: only vary params 2, 4, 6, and 8

# Define parameter ranges (min, max)
maxmin = [
    0.95  0.95;  # beta (Discount factor)
    0.10  0.20;  # delta (Depreciation rate)
    0.50  0.85;  # rho_e (Persistence of exchange rate shock)
    0.2   0.80;  # sigma_e (Volatility of exchange rate shock)
    0.40  0.80;  # nu (Share parameter for nondurable consumption)
    2.00  2.50;  # gamma (Risk aversion)
    0.05  0.80;  # f (Adjustment fixed cost)
    100   500;   # w (Wage)
    0.1   0.9;   # chi (Required maintenance)
    3000   5000;   # pd (Price of durables)
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

        answ = valfun(ppp);
        adjust_result=answ.adjust_result;
        noadjust_result=answ.noadjust_result;

        # Compute the moments using the updated parameters
        simdata = simmodel(answ)
        moments = makemoments(simdata, ppp, adjust_result, noadjust_result)

        # Store the parameter values and moments
        allparams[ivary, iparam] = glop
        allmoms[ivary, iparam, :] = moments

        # Compute adjustment gaps for the current parameter
        gap, f_x, x_values, h_x, I_d = adjustment_gaps(adjust_result, noadjust_result)
        gap_vec = vec(gap)

        # Call plotgaps_comp inside the loop
        plotgaps_comp(x_values, f_x, h_x, gap_vec, pname[iparam], glop)

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