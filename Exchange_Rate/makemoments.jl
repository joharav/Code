function makemoments(simdata::NamedTuple, pea::Vector{Float64})
    # Initialize the output moments vector
    outmoms = zeros(sz.nmom)
    
    # Constants from `pea`
    beta = pea[1]
    delta = pea[2]
    f = pea[7]
    w = pea[10]
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
        if current_d[i, j] == d_next[i, j]
            c[i, j] = w + current_e[i, j] * current_a[i, j] * (1 + rr) + current_e[i, j] * current_p[i, j] * (1 - f) * (1 - delta) * current_d[i, j] - a_next[i, j]
        else
            c[i, j] = w + current_e[i, j] * current_a[i, j] * (1 + rr) + current_e[i, j] * current_p[i, j] * (1 - delta) * current_d[i, j] - a_next[i, j] - current_p[i, j] * d_next[i, j]
        end
    end

    # Calculate investment rates and their differences
    skay = size(d, 1)
    invest = (d[4:skay-0, :] .- (1.0 - delta) .* d[3:skay-1, :]) ./ d[3:skay-1, :]

    # Moments calculations
    mu_i = mean(vec(invest))
    var_i = var(vec(invest))
    mu_a = mean(vec(a))
    var_a = var(vec(a))
    mu_c = mean(vec(c))
    var_c = var(vec(c))
    avg_spell_length = mean(diff(findall(current_d .!= d_next)))
    prob_change = mean(current_d .!= d_next)

    # # Calculate the average spell length between changes of the durable good `d`
    # spell_lengths = []
    # change_count = 0
    # total_periods = 0

    # Threads.@threads for j in 1:sz.nFirms
    #     last_value = d[1, j]
    #     spell_length = 0
    
    #     Threads.@threads for i in max(1, sz.burnin-2):min(sz.nYears, size(d, 1))
    #         total_periods += 1
    #         if d[i, j] == last_value
    #             spell_length += 1
    #         else
    #             push!(spell_lengths, spell_length)
    #             change_count += 1
    #             spell_length = 1
    #             last_value = d[i, j]
    #         end
    #     end
    #     push!(spell_lengths, spell_length)  # Final spell length
    # end
    # avg_spell_length = mean(spell_lengths)
    # prob_change = change_count / total_periods

    # Populate outmoms
    outmoms[1] = mu_i
    outmoms[2] = var_i
    outmoms[3] = mu_a
    outmoms[4] = var_a
    outmoms[5] = mu_c
    outmoms[6] = var_c
    outmoms[7] = avg_spell_length  
    outmoms[8] = prob_change       
    # Calculate kurtosis
    kurt_i = kurtosis(vec(current_a))
    kurt_a = kurtosis(vec(a))
    kurt_c = kurtosis(vec(c))

    # Calculate ratios
    ratio_d_income = mean(vec(current_e .* current_d) ./ vec(w .+ current_e .* current_a .* (1 + rr)))
    ratio_d_wealth = mean(vec(current_e .* current_d) ./ vec(current_e .* current_a .+ current_e .* current_d))
    ratio_d_consumption = mean(vec(current_e .* current_d) ./ vec(c))

    # Populate outmoms
    outmoms[1] = mu_i
    outmoms[2] = var_i
    outmoms[3] = mu_a
    outmoms[4] = var_a
    outmoms[5] = mu_c
    outmoms[6] = var_c
    outmoms[7] = avg_spell_length  
    outmoms[8] = prob_change
    outmoms[9] = kurt_i
    outmoms[10] = kurt_a
    outmoms[11] = kurt_c
    outmoms[12] = ratio_d_income
    outmoms[13] = ratio_d_wealth
    outmoms[14] = ratio_d_consumption

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
        println("Average spell length between changes in durables: $avg_spell_length\n")
        println("Probability of change in durables: $prob_change\n")
        println("Kurtosis of rate of investment: $kurt_i\n")
        println("Kurtosis of assets: $kurt_a\n")
        println("Kurtosis of nondurable consumption: $kurt_c\n")
        println("Ratio of durable holdings to income: $ratio_d_income\n")
        println("Ratio of durable holdings to wealth: $ratio_d_wealth\n")
        println("Ratio of durable holdings to consumption: $ratio_d_consumption\n")
        println("----------------------------------------------------------")
    end 

    return outmoms::Vector{Float64}
end