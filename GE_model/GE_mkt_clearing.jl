function GE_errors(pea)
    R_low = 1.0
    R_high = 1.0 / sz.beta

    for iter in 1:max_iter
        R = 0.5 * (R_low + R_high)
        K_demand = (sz.alpha / (R -1 + sz.delta))^(1 / (1 - sz.alpha))  # Cobb-Douglas demand
        WW = (1-sz.alpha)*K_demand^sz.alpha


        # Update parameters, solve HH block
        pea_ge = copy(pea)
        theta = pea[14]           # Share of savings in dollars
        R_star = pea[15]          # Foreign return
        E_e = median(grids.ex)    # Expected exchange rate

        R_eff = R * (1 - theta) + R_star * E_e * theta
        pea_ge[1] = 1 / R_eff  # beta = 1 / R
        pea_ge[8] = WW      # wage
        pd= pea_ge[10]  # price of durable goods

        # Solve household problem
        answ = valfun(pea_ge)

        # Simulate model
        simdata = simmodel(answ)

        # ----- Market 1: Asset market -----
        a_sim                   = simdata.a[sz.burnin-2:sz.nYears, :]
        A_supply                = mean(sum(a_sim, dims=2))  # household savings
        error_K                 = K_demand - A_supply

        if error_K <0
            R_low = R
        else
            R_high = R
        end

       # ----- Market 2: Durable market -----

       X = mean(sum(simdata.d, dims=2))                     # Durable demand and supply
       E_e = median(grids.ex)                               # Avg exchange rate


        # ----- Market 3: Labor market -----

        Lf = K_demand^sz.alpha
        Ld = (E_e * pd * X) / WW
        L_total = Lf + Ld
        
        # Calculate error in K market
        if abs(error_K) < tol
            return R, K_demand, WW, simdata
        end


    end
    println("Firm labor = $(Lf), Durable labor = $(Ld), Total labor demand = $(L_total), Return = $(R), K_demand = $(K_demand), Wage = $(WW)")

    error("Bisection failed to converge after $max_iter iterations")

end
