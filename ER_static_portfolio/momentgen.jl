# ==========================================================================
# 4D MODEL: Moment generation wrapper
# Main entry point: takes parameters, solves model, simulates, returns moments
# ==========================================================================

function momentgen(p::Vector{Float64})
    commence = time()
    
    # Solve model (both regimes)
    answ = valfun(p)
    
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
    moms = makemoments(simdata, p; per_year=4, shock=false)
    
    total_time = time() - commence
    if settings.verbose
        println("Total elapsed time: ", round(total_time, digits=1), " seconds")
    end
    
    # Optional: decision rules visualization
    if settings.compstat == false && settings.verbose == true
        decision_rules(answ)
    end
    
    # Optional: welfare analysis
    if settings.welfare
        run_welfare_analysis(p, answ, simdata)
    end
    
    return moms::Vector{Float64}
end


# ==========================================================================
# Helper for welfare analysis (placeholder - implement if needed)
# ==========================================================================
function run_welfare_analysis(p::Vector{Float64}, answ::NamedTuple, simdata::NamedTuple)
    println("\n=== Welfare Analysis ===")
    
    # Baseline welfare
    v_mean = mean(answ.v)
    println("Mean value function: ", round(v_mean, digits=4))
    
    # Counterfactual 1: Higher fixed cost
    pe_B = copy(p)
    pe_B[7] = 1.0  # f = 1
    # results_B = welfare_compare(p, pe_B)
    
    # Counterfactual 2: Fixed ER
    pe_C = copy(p)
    pe_C[4] = 0.0  # sigma_e = 0
    # results_C = welfare_compare(p, pe_C)
    
    println("(Detailed welfare analysis not yet implemented in 4D)")
end


# ==========================================================================
# Decision rules visualization (placeholder)
# ==========================================================================
function decision_rules(answ::NamedTuple)
    if !settings.verbose
        return
    end
    
    grids = answ.g
    pol = answ.pol
    
    println("\n=== Policy Function Summary ===")
    
    # Summary statistics at median states
    ie_mid = div(sz.ne, 2) + 1
    iy_mid = div(sz.ny, 2) + 1
    
    println("At median (e, y) state:")
    println("  Mean w': ", round(mean(pol.w[ie_mid, iy_mid, :, :]), digits=2))
    println("  Mean d': ", round(mean(pol.d[ie_mid, iy_mid, :, :]), digits=2))
    println("  Mean s:  ", round(mean(pol.s[ie_mid, iy_mid, :, :]), digits=4))
    
    # Dollar share statistics
    s_all = vec(pol.s)
    println("\nDollar share (s) distribution:")
    println("  Mean: ", round(mean(s_all), digits=4))
    println("  Std:  ", round(std(s_all), digits=4))
    println("  Min:  ", round(minimum(s_all), digits=4))
    println("  Max:  ", round(maximum(s_all), digits=4))
end
