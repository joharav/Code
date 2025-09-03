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
            #Large Fixed Cost on Durables
            pe_B = copy(p)
            pe_B[7] = 1
            pe_B[11] = 1
            results = welfare_full_summary(p, pe_B)


            #Large Fixed Cost on Durables
            pe_B2 = copy(p)
            pe_B2[7] = 1
            resultsB2 = welfare_full_summary(p, pe_B2)

            #Large Fixed Cost on Durables
            pe_B3 = copy(p)
            pe_B3[11] = 1
            resultsB3 = welfare_full_summary(p, pe_B3)

            #Fixed ER
            pe_C = copy(p)
            pe_C[4] = 0
            resultsC = welfare_full_summary(p, pe_C)

            #High Volatility
            pe_D = copy(p)
            pe_D[4] = 2
            resultsD = welfare_full_summary(p, pe_D)


            # Collect all results
            cases = ["Durable FCost 1", "Durable FCost 2", "Durable FCost 3", "Fixed ER", "High Volatility"]
            cev_BA_vals = [results.cev_BA, resultsB2.cev_BA, resultsB3.cev_BA, resultsC.cev_BA, resultsD.cev_BA]
            cev_AB_vals = [results.cev_AB, resultsB2.cev_AB, resultsB3.cev_AB, resultsC.cev_AB, resultsD.cev_AB]
            lambda_vals = [results.λ_composite, resultsB2.λ_composite, resultsB3.λ_composite, resultsC.λ_composite, resultsD.λ_composite]
            λ_c_only_vals = [results.λ_c_only, resultsB2.λ_c_only, resultsB3.λ_c_only, resultsC.λ_c_only, resultsD.λ_c_only]
            keepDistAB_vals = [results.keepDistAB, resultsB2.keepDistAB, resultsB3.keepDistAB, resultsC.keepDistAB, resultsD.keepDistAB]
            keepDistBA_vals = [results.keepDistBA, resultsB2.keepDistBA, resultsB3.keepDistBA, resultsC.keepDistBA, resultsD.keepDistBA]
            acrossSS_vals = [results.acrossSS, resultsB2.acrossSS, resultsB3.acrossSS, resultsC.acrossSS, resultsD.acrossSS]
            
            # Create a DataFrame
            df = DataFrame(
                Case = cases,
                CEV_BA = cev_BA_vals,
                #CEV_AB = cev_AB_vals,
                Lambda = lambda_vals,
                #Lambda_c = λ_c_only_vals,    
                WelfareChange_A_to_B_keep_dist = keepDistAB_vals,
              #  WelfareChange_B_to_A_keep_dist = keepDistBA_vals,
                WelfareChange_across_SS = acrossSS_vals
            )
            
            # Export to CSV
            CSV.write("Output/Welfare_Comparison.csv", df)
            

        end


    else
        moms = -100.0*ones(sz.nmom)
    end
    
    return moms::Vector{Float64}
end