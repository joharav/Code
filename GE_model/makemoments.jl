function makemoments(simdata::NamedTuple, pea::Vector{Float64}; shock::Bool = false)
    # Initialize the output moments vector
    outmoms = zeros(sz.nmom)
    
    # Constants from `pea`
    beta  = pea[1]
    w     = pea[8]
    pd    = pea[10]
    theta = pea[16]
    R_star= pea[17]
    R     = 1/beta
    R_eff = (1 - theta) * R + theta * R_star

    # Extract variables from simulation data
    a                   = simdata.a[sz.burnin-2:sz.nYears, :]
    a_state             = simdata.a[sz.burnin-3:sz.nYears-1, :]
    d                   = simdata.d[sz.burnin-2:sz.nYears, :]
    d_state             = simdata.d[sz.burnin-3:sz.nYears-1, :]
    ex                  = simdata.ex[sz.burnin-2:sz.nYears, :]
    c                   = simdata.c[sz.burnin-2:sz.nYears, :]
    d_adjust            = simdata.d_adjust[sz.burnin-2:sz.nYears, :]
    adjust_indicator    = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]
    zz                  = simdata.zz[sz.burnin-2:sz.nYears, :]

    # Calculate the gaps
    adjustment_indicator = vec(adjust_indicator)
    gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio =adjustment_gaps_sim(d_state,d_adjust,adjustment_indicator)
    d_invest = 100*(d .- d_state)./d_state
    a_change = 100*(a .- a_state)./a_state
    c_change = 100 * (c[2:end] .- c[1:end-1]) ./ c[1:end-1]
    mean_adjust_size = mean(vec(d[adjust_indicator .== 1]) .- vec(d_state[adjust_indicator .== 1]))


    # Moments calculations
    mu_d = mean(vec(d_invest))
    var_d = var(vec(d_invest))
    mu_a = mean(vec(a_change))
    var_a = var(vec(a_change))
    mu_c = mean(vec(c_change))
    var_c = var(vec(c_change))
    mu_d1 = mean(vec(d_state))
    var_d1 = var(vec(d_state))

    # Calculate ratios
    ratio_d_income = (vec(pd.* ex .* d) ./ vec(w .* zz .+ a_state .* R_eff ))
    ratio_d_wealth = (vec(pd.*ex .* d) ./ vec(a_state .* R_eff .+ pd*ex .* d_state))

    # Calculate the ratio
    ratio_d_consumption = (vec(pd.* ex .* d) ./ vec(c))
   
    mu_d_income = mean(ratio_d_income)
    mu_d_wealth = mean(ratio_d_wealth) 
    mu_d_c      = mean(ratio_d_consumption)

    # Calculate distributions using KDE
    kde_ratio_d_income          = kde(ratio_d_income)
    kde_ratio_d_wealth          = kde(ratio_d_wealth)
    kde_ratio_d_consumption     = kde(ratio_d_consumption)

    # Distribution of simulated KDE 
    f_d_income      = kde_ratio_d_income.density
    f_d_wealth      = kde_ratio_d_wealth.density
    f_d_consumption = kde_ratio_d_consumption.density

    x_values_d_income       = collect(kde_ratio_d_income.x)
    x_values_d_wealth       = collect(kde_ratio_d_wealth.x)
    x_values_d_consumption  = collect(kde_ratio_d_consumption.x)
    
    #Dispersion measures
    disp_d_income           =compute_dispersion(ratio_d_income)
    disp_d_wealth           =compute_dispersion(ratio_d_wealth)
    disp_d_c                =compute_dispersion(ratio_d_consumption)
    disp_d                  =compute_dispersion(vec(d))

    #Interquartile range
    IQR_d_income            =disp_d_income[2]
    IQR_d_wealth            =disp_d_wealth[2]
    IQR_d_c                 =disp_d_c[2]
    IQR_d                   =disp_d[2]

    #90th to 10th percentile ratio
    p90_10_d_income         =disp_d_income[3]
    p90_10_d_wealth         =disp_d_wealth[3]
    p90_10_d_c              =disp_d_c[3]
    p90_10_d                =disp_d[3]



    # Populate outmoms
    outmoms[1]  = mu_d
    outmoms[2]  = var_d
    outmoms[3]  = mu_a
    outmoms[4]  = var_a
    outmoms[5]  = mu_c
    outmoms[6]  = var_c
    outmoms[7]  = mu_d_income
    outmoms[8]  = mu_d_wealth
    outmoms[9]  = mu_d_c
    outmoms[10] = mu_gap
    outmoms[11] = var_gap
    outmoms[12] = mean_adjust_size
    outmoms[13] = adjustment_ratio
    if settings.compstat==false
        plotgaps(x_values, f_x, h_x, gap_vec; shock=shock)
        plotdensities(x_values_d_income, f_d_income, "f_income"; shock=shock)
        plotdensities(x_values_d_wealth, f_d_wealth, "f_wealth"; shock=shock)
        plotdensities(x_values_d_consumption, f_d_consumption, "d_c"; shock=shock)
        plot_aggregates(simdata)
        d_adjust_time_size(simdata)
    end

    println("Adjustment Ratio: $adjustment_ratio\n")    


    if settings.verbose && settings.irfsshock==false

        println("----------------------------------------------------------")
        println("\nStatistics:\n")
        println("Average durables: $mu_d\n")
        println("Variance of the durable holdings: $var_d\n")
        println("Average assets: $mu_a\n")
        println("Variance of assets: $var_a\n")
        println("Average nondurable consumption: $mu_c\n")
        println("Variance of nondurable consumption: $var_c\n")
        println("Ratio of durable holdings to income: $mu_d_income\n")
        println("Ratio of durable holdings to wealth: $mu_d_wealth\n")
        println("Ratio of durable holdings to consumption: $mu_d_c\n")
        println("Average gap: $mu_gap\n")
        println("Variance of gap: $var_gap\n")
        println("Mean adjustment size: $mean_adjust_size\n")
        println("Adjustment Ratio: $adjustment_ratio\n")    
        println("----------------------------------------------------------")
    
        println("----------------------------------------------------------")
        println("\nInterquartile ratio:\n")
        println("Durable stock to income: $IQR_d_income\n")
        println("Durable stock to wealth: $IQR_d_wealth\n")
        println("Durable to consumption: $IQR_d_c\n")
        println("Durables: $IQR_d\n")
        println("\nPercentile 90th to 10th ratio:\n")
        println("Durable stock to income: $p90_10_d_income\n")
        println("Durable stock to wealth: $p90_10_d_wealth\n")
        println("Durable to consumption: $p90_10_d_c\n")
        println("Durables: $p90_10_d\n")
        println("----------------------------------------------------------")

    end

        return outmoms::Vector{Float64}, x_values, f_x, h_x

end