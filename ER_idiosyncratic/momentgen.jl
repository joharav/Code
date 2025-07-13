function momentgen(p::Vector{Float64})
    # ============ Run stuff ===================================
    commence = time()
    answ = valfun(p)
    arret = time()
    println("elapse of time in seconds = ", arret-commence)

    if answ.e == 0
        simdata = simmodel(answ)
        ergodic_dist = compute_ergodic(answ)
        pe_baseline = copy(p)
        pe_baseline[16] = 0.0
        result0 = valfun(pe_baseline)
        v_base = vec(mean(result0.v, dims=(1,2,3)))
        v_current = vec(mean(answ.v, dims=(1,2,3)))
        cev = compute_cev(v_base, v_current, p)

        # ============ MOMENTS ===================================
        moms, x_values, f_x, h_x, _, _ = makemoments(simdata, p; shock=false,cev=cev)
        if settings.compstat==false
            decision_rules(answ)
        end        
        if settings.irfsshock
            # Simulate shock 
            simdata_irf = simmodel_girf(answ, Int(sz.nYears/2))
            
            # Get moments for simulation with shock
            moms_shock, x_values_shock, f_x_shock, h_x_shock,  _, _ = makemoments(simdata_irf, p; shock=true,cev=0)
            
            # Create GIRF plots
            girf = girf_plots(simdata_irf, simdata)
            cirf_c = compute_cirf(vec(girf[1]), 8, "c")
            cirf_d = compute_cirf(vec(girf[2]), 8, "d")
            cirf_a = compute_cirf(vec(girf[3]), 8, "a")
            plotgaps_shock(x_values, f_x, h_x, x_values_shock, f_x_shock, h_x_shock)
        end
        # ============ WELFARE COMPARISON ==============================
        if settings.welfare
            println("Computing welfare comparison...")
            run_batch()
            run_all_batches()
        end


    else
        moms = -100.0*ones(sz.nmom)
    end
    
    return moms::Vector{Float64}
end