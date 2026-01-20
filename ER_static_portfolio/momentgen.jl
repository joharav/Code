# ==========================================================================
# 4D MODEL: Moment generation wrapper
# Main entry point: takes parameters, solves model, simulates, returns moments
# ==========================================================================

function momentgen(p::Vector{Float64}; grid_builder = makegrids)
    commence = time()
    
    # Solve model (both regimes)
    answ = valfun(p; grid_builder = makegrids)
    
    solve_time = time()
    if settings.verbose
        println("\nSolve time: ", round(solve_time - commence, digits=1), " seconds")
    end
    
    # Simulate panel
    simdata = simmodel(answ)
    
    sim_time = time()
    if settings.verbose
        println("Simulation time: ", round(sim_time - solve_time, digits=1), " seconds")
    end
    
    # Compute moments
    moms, x_values, f_x, h_x, gap_vec = makemoments(simdata, p; per_year=4, shock=false)
    
    total_time = time() - commence
    if settings.verbose
        println("Total elapsed time: ", round(total_time, digits=1), " seconds")
    end
    
    # Optional: decision rules visualization
    if settings.compstat == false && settings.verbose == true
        decision_rules(answ)
        plot_aggregates(simdata)
        plot_simulated_d_and_a_by_state(simdata)
        d_adjust_time_size(simdata, p)
        compute_d_star_panel(simdata, answ)
        plot_distribution_panels(simdata, p)
        printstuff(answ)
        plotgaps(x_values::Vector{Float64}, f_x::Vector{Float64}, h_x::Vector{Float64}, gap_vec::Vector{Float64})
    end
    

    if settings.irfsshock
        # Simulate shock 
        simdata_irf = simmodel_girf(answ, Int(sz.nYears/2))
        
        # Get moments for simulation with shock
        moms_shock, x_values_shock, f_x_shock, h_x_shock,gap_vec_shock = makemoments(simdata_irf, p; shock=true)
        
        # Create GIRF plots
        girf = girf_plots(simdata_irf, simdata)
        cirf_c = compute_cirf(vec(girf[1]), 8; name="c")
        cirf_d = compute_cirf(vec(girf[2]), 8; name="d")
        cirf_a = compute_cirf(vec(girf[3]), 8; name="a")
        plotgaps_shock(x_values, f_x, h_x, x_values_shock, f_x_shock, h_x_shock)
    end


    # Optional: welfare analysis
    if settings.welfare
        ergodic_dist = compute_ergodic(answ)
        run_welfare_analysis(p, answ)
    end
    
    return moms::Vector{Float64}
end



# ==========================================================================
# Helper for welfare analysis (placeholder - implement if needed)
# ==========================================================================
function run_welfare_analysis(p::Vector{Float64}, answ::NamedTuple)
    println("\n=== Welfare Analysis ===")
    
    # Baseline welfare
    v_mean = mean(answ.v)
    println("Mean value function: ", round(v_mean, digits=4))
    
    # Counterfactual 1: Higher fixed cost
    pe_B = copy(p)
    pe_B[7] = 1.0  # f = 1
    pe_B[17] = 1
    results_B = welfare_summary(p, pe_B)


    #Large Fixed Cost on Durables
    pe_B2 = copy(p)
    pe_B2[7] = 1
    resultsB2 = welfare_summary(p, pe_B2)

    #Large Fixed Cost on Durables
    pe_B3 = copy(p)
    pe_B3[17] = 1
    resultsB3 = welfare_summary(p, pe_B3)

    
    # Counterfactual 2: Fixed ER
    pe_C = copy(p)
    pe_C[4] = 0.0001  # sigma_e = 0
    resultsC = welfare_summary(p, pe_C)

    # Counterfactual 3: High Volatility
    pe_D = copy(p)
    pe_D[4] = 2
    resultsD = welfare_summary(p, pe_D)

    # Counterfactual 4: High Depreciation
    pe_E = copy(p)
    pe_E[2] = 0.15
    resultsE = welfare_summary(p, pe_E)

    # Counterfactual 5: High Kappa
    pe_F = copy(p)
    pe_F[11] = 2
    resultsF = welfare_summary(p, pe_F)

    # Collect all results
    cases = ["Durable FCost 1", "Durable FCost 2", "Durable FCost 3", "Fixed ER", "High Volatility"]
    cev_BA_vals = [results.cev_BA, resultsB2.cev_BA, resultsB3.cev_BA, resultsC.cev_BA, resultsD.cev_BA,resultsE.cev_BA,resultsF.cev_BA]
    keepDistAB_vals = [results.keepDistAB, resultsB2.keepDistAB, resultsB3.keepDistAB, resultsC.keepDistAB, resultsD.keepDistAB, resultsE.keepDistAB, resultsF.keepDistAB]
    acrossSS_vals = [results.acrossSS, resultsB2.acrossSS, resultsB3.acrossSS, resultsC.acrossSS, resultsD.acrossSS, resultsE.acrossSS, resultsF.acrossSS]
    
    # Create a DataFrame
    df = DataFrame(
        Case = cases,
        CEV_BA = cev_BA_vals,
        WelfareChange_A_to_B_keep_dist = keepDistAB_vals,
        WelfareChange_across_SS = acrossSS_vals
    )
    
    # Export to CSV
    CSV.write("Output/Welfare_Comparison.csv", df)
    
end

