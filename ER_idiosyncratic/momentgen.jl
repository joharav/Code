function momentgen(p::Vector{Float64})
    # ============ Run stuff ===================================
    commence = time()
    answ = valfun(p)
    arret = time()
    println("elapse of time in seconds = ", arret-commence)

    if answ.e == 0
        simdata = simmodel(answ)
        moms = makemoments(simdata, p)
        
        if settings.irfsshock
            # Simulate shock 
            simdata_irf = simmodel_girf(answ, Int(sz.nYears/2))
            
            # Get moments for simulation with shock
            moms_shock = makemoments(simdata_irf, p)
            
            # Create GIRF plots
            girf = girf_plots(simdata_irf, simdata)
            cirf_c = compute_cirf(vec(girf[1]), 8, "c")
            cirf_d = compute_cirf(vec(girf[2]), 8, "d")
            cirf_a = compute_cirf(vec(girf[3]), 8, "a")
            
        end
    else
        moms = -100.0*ones(sz.nmom)
    end
    
    return moms::Vector{Float64}
end