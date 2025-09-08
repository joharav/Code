function momentgen(p::Vector{Float64})
    # ============ Run stuff ===================================
    commence = time()
    answ = valfun(p)
    arret = time()
    println("elapse of time in seconds = ", arret-commence)

    if answ.e == 0
        simdata = simmodel(answ)
        ergodic_dist = compute_ergodic(answ)


        # ============ MOMENTS ===================================
        moms, x_values, f_x, h_x, _, _ = makemoments(simdata, p; shock=false)
        plot_distribution_panels(simdata, p)
        if settings.compstat==false && settings.verbose==true
            decision_rules(answ)
        end        
        if settings.irfsshock
            # Simulate shock 
            simdata_irf = simmodel_girf(answ, Int(sz.nYears/2))
            
            # Get moments for simulation with shock
            moms_shock, x_values_shock, f_x_shock, h_x_shock,  _, _ = makemoments(simdata_irf, p; shock=true)
            
            # Create GIRF plots
            girf = girf_plots(simdata_irf, simdata)
            cirf_c = compute_cirf(vec(girf[1]), 8, "c")
            cirf_d = compute_cirf(vec(girf[2]), 8, "d")
            cirf_a = compute_cirf(vec(girf[3]), 8, "a")
            plotgaps_shock(x_values, f_x, h_x, x_values_shock, f_x_shock, h_x_shock)
        end
        # ============ WELFARE COMPARISON ==============================
        if settings.welfare
            #No Dollars
            pe_nd = copy(p)
            pe_nd[16] = 0.1
            results_nd = welfare_full_summary(p, pe_nd)



            #Large Fixed Cost on Durables
            pe_B = copy(p)
            pe_B[7] = 1
            pe_B[11] = 1
           # pe_B[16] = 0
            results = welfare_full_summary(p, pe_B)

            #No Durables
            pe_zero = copy(p)
            pe_zero[2] = 0.9
            results_zero = welfare_full_summary(p, pe_zero)

            #Fixed ER
            pe_C = copy(p)
            pe_C[4] = 0
          #  pe_C[16] = 0
            resultsC = welfare_full_summary(p, pe_C)

            #High Volatility
            pe_D = copy(p)
            pe_D[4] = 1.8
           # pe_D[16] = 0
            resultsD = welfare_full_summary(p, pe_D)


            # Collect all results
            cases = [
                "No dollar savings",
                "High adj. cost",
                "No durables",
                "Fixed ER",
                "High ER volatility"
            ]
            cev_BA_vals = [results_nd.cev_BA, results.cev_BA, results_zero.cev_BA ,resultsC.cev_BA, resultsD.cev_BA]
          #  lambda_vals = [results_nd.λ_composite, results.λ_composite, results_zero.λ_composite, resultsC.λ_composite, resultsD.λ_composite]
           # keepDistAB_vals = [results_nd.keepDistAB, results.keepDistAB, results_zero.keepDistAB, resultsC.keepDistAB, resultsD.keepDistAB]
            #acrossSS_vals = [results_nd.acrossSS, results.acrossSS, results_zero.acrossSS, resultsC.acrossSS, resultsD.acrossSS]
            vol_c_vals = [results_nd.vol_c_with, results.vol_c_with, results_zero.vol_c_with, resultsC.vol_c_with, resultsD.vol_c_with]
            vol_x_vals = [results_nd.vol_x_with, results.vol_x_with, results_zero.vol_x_with, resultsC.vol_x_with, resultsD.vol_x_with]
           # vol_cno_vals = [results_nd.vol_c_no, results.vol_c_no, results_zero.vol_c_no, resultsC.vol_c_no, resultsD.vol_c_no]
           # vol_xno_vals = [results_nd.vol_x_no, results.vol_x_no, results_zero.vol_x_no, resultsC.vol_x_no, resultsD.vol_x_no]
                   
            # Create a DataFrame
            df = DataFrame(
                Case = cases,
                CEV_BA = cev_BA_vals,
              #  Lambda = lambda_vals,
               # WelfareChange_A_to_B_keep_dist = keepDistAB_vals,
               # WelfareChange_across_SS = acrossSS_vals,
                ConsumptionVolatility = vol_c_vals,
                DurVolatility = vol_x_vals
               # ConsumptionVolatility_no_dollars = vol_cno_vals,
                #DurVolatility_no_dollars = vol_xno_vals
            )
            
            # Export to CSV
            CSV.write("Output/Welfare_Comparison.csv", df)
            

        end


    else
        moms = -100.0*ones(sz.nmom)
    end
    
    return moms::Vector{Float64}
end