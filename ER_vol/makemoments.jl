function makemoments(simdata::NamedTuple, pea::Vector{Float64})
    # Initialize the output moments vector
    outmoms = zeros(sz.nmom)
    
    # Constants from `pea`
    beta  = pea[1]
    delta = pea[2]
    f     = pea[5]
    w     = pea[8]
    chi   = pea[9]
    pd    = pea[10]
    rr = 1/beta -1

    # Extract variables from simulation data
    a = simdata.a[sz.burnin-2:sz.nYears, :]
    d = simdata.d[sz.burnin-2:sz.nYears, :]
    ex = simdata.ex[sz.burnin-2:sz.nYears, :]
    d_adjust = simdata.d_adjust[sz.burnin-2:sz.nYears, :]
    adjust_indicator = simdata.adjust_indicator[sz.burnin-2:sz.nYears, :]

    # Prepare slices of next period's asset and durable values
    a_next = a[2:end, :]
    d_next = d[2:end, :]

    # Slice current period's variables correctly
    current_a = a[1:end-1, :]
    current_d = d[1:end-1, :]
    current_e = ex[1:end-1, :]
    d_a = d_adjust[1:end-1, :]
    adjust_indicator= adjust_indicator[1:end-1, :]
    # Calculate consumption
    # Initialize `c` with the same size as `current_a`
    c = zeros(size(current_a))
    # Compute consumption using a loop
    for i in 1:sz.nYears-(sz.burnin-2), j in 1:sz.nFirms
        if d_next[i, j] == current_d[i, j] * (1 - delta * (1 - chi))  
            c[i, j] = w + current_e[i, j] * current_a[i, j] * (1 + rr) - current_e[i, j] * pd * delta * chi * current_d[i, j] - current_e[i, j] * a_next[i, j]
        else
            c[i, j] = w + current_e[i, j] * current_a[i, j] * (1 + rr) + current_e[i, j] * pd * (1 - delta) * (1-f) * current_d[i, j] 
            - current_e[i, j] * a_next[i, j] - current_e[i, j] * pd * d_next[i, j]
        end
    end
    adjustment_indicator = vec(adjust_indicator)
    gap_vec, f_x, x_values, h_x, I_d, mu_gap, var_gap, adjustment_ratio =adjustment_gaps_sim(current_d,d_a,adjustment_indicator)

    # Moments calculations
    mu_d = mean(vec(current_d))
    var_d = var(vec(current_d))
    mu_a = mean(vec(current_a))
    var_a = var(vec(current_a))
    mu_c = mean(vec(c))
    var_c = var(vec(c))

    # Calculate ratios
    ratio_d_income = (vec(pd.* current_e .* current_d) ./ vec(w .+ current_e .* current_a .* (1 + rr) ))
    ratio_d_wealth = (vec(pd.*current_e .* current_d) ./ vec(current_e .* current_a .* (1 + rr) .+ pd*current_e .* current_d))
   
    # Ensure the vectors have the same length
    min_length = min(length(vec(pd.*current_e .* current_d)), length(vec(c)))

    # Truncate the vectors to the minimum length
    truncated_numerator = vec(pd.*current_e .* current_d)[1:min_length]
    truncated_denominator = vec(c)[1:min_length]

    # Calculate the ratio
    ratio_d_consumption = (truncated_numerator ./ truncated_denominator)
   
    mu_d_income = mean(ratio_d_income)
    mu_d_wealth = mean(ratio_d_wealth) 
    mu_d_c      = mean(ratio_d_consumption)

    # Calculate distributions using KDE
    kde_ratio_d_income = kde(ratio_d_income)
    kde_ratio_d_wealth = kde(ratio_d_wealth)
    kde_ratio_d_consumption = kde(ratio_d_consumption)

    # Distribution of simulated KDE 
    f_d_income      = kde_ratio_d_income.density
    f_d_wealth      = kde_ratio_d_wealth.density
    f_d_consumption = kde_ratio_d_consumption.density

    x_values_d_income = collect(kde_ratio_d_income.x)
    x_values_d_wealth = collect(kde_ratio_d_wealth.x)
    x_values_d_consumption = collect(kde_ratio_d_consumption.x)

    # Populate outmoms
    outmoms[1] = mu_d
    outmoms[2] = var_d
    outmoms[3] = mu_a
    outmoms[4] = var_a
    outmoms[5] = mu_c
    outmoms[6] = var_c
    outmoms[7] = mu_d_income
    outmoms[8] = mu_d_wealth
    outmoms[9] = mu_d_c
    outmoms[10] = mu_gap
    outmoms[11] = var_gap
    outmoms[12] = I_d
    outmoms[13] = adjustment_ratio


    if settings.verbose 

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
        println("Aggregate durable expenditures: $I_d\n")
        println("Adjustment Ratio: $adjustment_ratio\n")    
        println("----------------------------------------------------------")


        plotgaps(x_values, f_x, h_x, gap_vec)
        # Plot the distribution f(x) separately
        plotdensities(x_values_d_income, f_d_income, "f_income")
        plotdensities(x_values_d_wealth, f_d_wealth, "f_wealth")
        plotdensities(x_values_d_consumption, f_d_consumption, "d_c")

    end

    # Optionally print statistics
    if settings.compstat
        return outmoms::Vector{Float64}, x_values, f_x, h_x, gap_vec
    else 
        return outmoms::Vector{Float64}
    end 
end