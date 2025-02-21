function makemoments(simdata::NamedTuple, pea::Vector{Float64}, adjust_result::NamedTuple, noadjust_result::NamedTuple)
    # Initialize the output moments vector
    outmoms = zeros(sz.nmom)
    
    # Constants from `pea`
    beta  = pea[1]
    delta = pea[2]
    f     = pea[7]
    w     = pea[10]
    chi   = pea[11]
    rr = 1/beta -1

    # Extract variables from simulation data
    a = simdata.a[sz.burnin-2:sz.nYears, :]
    d = simdata.d[sz.burnin-2:sz.nYears, :]
    p = simdata.p[sz.burnin-2:sz.nYears, :]
    ex = simdata.ex[sz.burnin-2:sz.nYears, :]

    # Prepare slices of next period's asset and durable values
    a_next = a[2:end, :]
    d_next = d[2:end, :]

    # Slice current period's variables correctly
    current_a = a[1:end-1, :]
    current_d = d[1:end-1, :]
    current_p = p[1:end-1, :]
    current_e = ex[1:end-1, :]


    # Calculate consumption
    # Initialize `c` with the same size as `current_a`
    c = zeros(size(current_a))

    # Compute consumption using a loop
    for i in 1:sz.nYears-sz.burnin-2, j in 1:sz.nFirms
        if d_next[i, j] == current_d[i, j] * (1 - delta * (1 - chi))  
            c[i, j] = w + current_e[i, j] * current_a[i, j] * (1 + rr) - current_e[i, j] * current_p[i, j] * delta * chi * current_d[i, j] - current_e[i, j] * a_next[i, j]
        else
            c[i, j] = w + current_e[i, j] * current_a[i, j] * (1 + rr) + current_e[i, j] * current_p[i, j] * (1 - delta) * (1-f) * current_d[i, j] - current_e[i, j] * a_next[i, j] - current_e[i, j] * current_p[i, j] * d_next[i, j]
        end
    end

    # Calculate investment rates and their differences
    skay = size(d, 1)
    invest = (d[4:skay-0, :] .- (1.0 - delta) .* d[3:skay-1, :]) ./ d[3:skay-1, :]

    # Moments calculations
    mu_i = mean(vec(invest))
    var_i = var(vec(invest))
    mu_a = mean(vec(current_a))
    var_a = var(vec(current_a))
    mu_c = mean(vec(c))
    var_c = var(vec(c))

    # Calculate ratios
    ratio_d_income = mean(vec(current_e .* current_d) ./ vec(w .+ current_e .* current_a .* (1 + rr)))
    ratio_d_wealth = mean(vec(current_e .* current_d) ./ vec(current_e .* current_a .* (1 + rr) .+ current_e .* current_d))
   
    # Ensure the vectors have the same length
    min_length = min(length(vec(current_e .* current_d)), length(vec(c)))

    # Truncate the vectors to the minimum length
    truncated_numerator = vec(current_e .* current_d)[1:min_length]
    truncated_denominator = vec(c)[1:min_length]

    # Calculate the ratio
    ratio_d_consumption = mean(truncated_numerator ./ truncated_denominator)
   
    # Adjustment Gaps and distributions

    gap, f_x, x_values, h_x, I_d = adjustment_gaps(adjust_result,noadjust_result)
    gap_vec=vec(gap)
    mu_gap = mean(gap)
    var_gap = var(gap)
    mu_hx = mean(h_x)
    var_hx = var(h_x)
    I_d=I_d


    # Populate outmoms
    outmoms[1] = mu_i
    outmoms[2] = var_i
    outmoms[3] = mu_a
    outmoms[4] = var_a
    outmoms[5] = mu_c
    outmoms[6] = var_c
    outmoms[7] = ratio_d_income
    outmoms[8] = ratio_d_wealth
    outmoms[9] = ratio_d_consumption
    outmoms[10] = mu_gap
    outmoms[11] = var_gap
    outmoms[12] = mu_hx
    outmoms[13] = var_hx
    outmoms[14] = I_d


    # Optionally print statistics
    if settings.compstat
        println("----------------------------------------------------------")
        println("\nStatistics:\n")
        println("Average rate of durable investment: $mu_i\n")
        println("Variance of the rate of durable investment: $var_i\n")
        println("Average assets: $mu_a\n")
        println("Variance of assets: $var_a\n")
        println("Average nondurable consumption: $mu_c\n")
        println("Variance of nondurable consumption: $var_c\n")
        println("Ratio of durable holdings to income: $ratio_d_income\n")
        println("Ratio of durable holdings to wealth: $ratio_d_wealth\n")
        println("Ratio of durable holdings to consumption: $ratio_d_consumption\n")
        println("Average gap: $mu_gap\n")
        println("Variance of gap: $var_gap\n")
        println("Average adjustment hazard: $mu_hx\n")
        println("Variance of adjustment hazard: $var_hx\n")
        println("Aggregate durable expenditures: $I_d\n")
        println("----------------------------------------------------------")
        plotgaps(x_values, f_x, h_x,gap_vec)
    end 

    return outmoms::Vector{Float64}
end