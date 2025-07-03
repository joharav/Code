function GE_errors(pea)
    beta= pea[1]  # Discount factor
    delta = pea[2]  # Depreciation rate
    alpha = pea[18]  # Capital share
    grids = makegrids(pea)
    R_low = 0.5   # Lower bound for interest rate
    R_high = 5

    for iter in 1:sz.maxiter
        R = 0.5 * (R_low + R_high)
        K_demand = (alpha / (R -1 + delta))^(1 / (1 - alpha))  # Cobb-Douglas demand
        WW = (1-alpha)*K_demand^alpha

        E_e = median(grids.ex)    # Expected exchange rate

        # Update parameters, solve HH block
        pea_ge = copy(pea)
        theta = pea[16]           # Share of savings in dollars
        R_star = pea[17]          # Foreign return
        R_eff = R * (1 - theta) + R_star * theta

        pea_ge[1] = 1 / R_eff  # beta = 1 / R
        pea_ge[8] = WW      # wage
        pd= pea_ge[10]  # price of durable goods

        println("GE Iteration $iter: R = $R, K_demand = $K_demand, W = $WW, R_eff = $R_eff")

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

        Lf = K_demand^alpha
        Ld = (E_e * pd * X) / WW
        L_total = Lf + Ld
        
       # Convergence check
    if abs(error_K) < sz.toler
        return (
            R = R,
            K = K_demand,
            W = W,
            R_eff = R_eff,
            E_e = E_e,
            simdata = simdata,
            A_supply = A_supply,
            X = X,
            Lf = Lf,
            Ld = Ld,
            L_total = L_total,
            converged = true
        )
        end

        # Update bisection bounds
        if error_K < 0
            R_low = R
        else
            R_high = R
        end
    end

    println("Bisection failed. Final stats:")
    println("Firm labor = $Lf, Durable labor = $Ld, Total = $L_total")
    println("R = $R, K = $K_demand, Wage = $W")

    return (converged = false)

    end
