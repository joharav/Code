function momentgen(p::Vector{Float64})
    start_time = time()
    


        if settings.GE

            # ============ GENERAL EQUILIBRIUM BLOCK ====================
            println("Solving GE equilibrium...")
            ge = GE_errors(p)

            if !ge.converged
                println("GE did not converge.")
                return -100.0 * ones(sz.nmom), nothing
            end

            # Update parameter vector with GE prices
            p_ge = copy(p)
            p_ge[1] = 1 / ge.R_eff   # Update β = 1 / R_eff
            p_ge[8] = ge.W           # Update wage
            p_ge[10] = p[10]         # Durable price stays? Or from GE?

            # Re-solve household problem
            answ = valfun(p_ge)
            simdata = ge.simdata  # Already simulated inside GE_errors

        else 
            # ============ PARTIAL EQUILIBRIUM BLOCK ====================

            println("Solving PE equilibrium...")
            answ = valfun(p)
            
            if answ.e != 0
                println("Value function iteration did not converge.")
                return -100.0 * ones(sz.nmom), answ
            end

            simdata = simmodel(answ)
            if settings.verbose
                println("PE simulation completed.")
            end

        end

    # ============ MOMENT GENERATION ===============================
    moms, x_values, f_x, h_x = makemoments(simdata, p; shock=false)

    if settings.compstat == false
        decision_rules(answ)
    end

    # ============ SHOCK SIMULATION (IRFs) =========================
    if settings.irfsshock
        simdata_irf = simmodel_girf(answ, Int(sz.nYears/2))
        moms_shock, x_vals_shock, f_x_shock, h_x_shock = makemoments(simdata_irf, p; shock=true)
        girf = girf_plots(simdata_irf, simdata)
        _ = compute_cirf(vec(girf[1]), 8, "c")
        _ = compute_cirf(vec(girf[2]), 8, "d")
        _ = compute_cirf(vec(girf[3]), 8, "a")
        plotgaps_shock(x_values, f_x, h_x, x_vals_shock, f_x_shock, h_x_shock)
    end

    # ============ WELFARE COMPARISON ==============================
    if settings.welfare
        println("Computing welfare comparison...")
        run_welfare_analysis()
    end

    # ============ MIT SHOCK (OPTIONAL) ============================
    if settings.mitshock
        println("Computing MIT shock dynamics...")
        mit_shock_dynamics(simdata)
    end

    elapsed = time() - start_time
    println("Elapsed time: ", round(elapsed, digits=2), " seconds.")

    return moms::Vector{Float64}, answ
end
