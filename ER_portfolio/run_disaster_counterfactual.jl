function run_disaster_counterfactual()
    # 1. Get baseline parameters
    p_base = baseline_params()
    
    # 2. Generate Disaster Grids
    # Example: 2.5% annual prob, 40% devaluation (log 1.4)
    g_disaster = makegrids_disaster(p_base; pi_annual=0.025, kappa_e=log(1.4))

    # 3. Solve Model with NEW grids
    # Pass the explicit grids to valfun (you might need to modify valfun to accept grids)
    # If your valfun calls makegrids internally, you need to create a version 
    # that accepts grids as an argument: valfun(params, grids)
    
    # Assuming you have/make a version: valfun(p, precomputed_grids)
    answ_disaster = valfun(p_base, g_disaster) 

    # 4. Simulate
    sim_disaster  = simmodel(answ_disaster)
    
    # 5. Analyze
    mech_disaster = mechanism_stats(sim_disaster, p_base)
    println("Dollar Share with Disaster Risk: ", mech_disaster.usd_share_mean)
end