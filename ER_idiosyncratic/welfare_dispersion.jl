function batch_welfare_and_dispersion(pe_base::Vector{Float64}, theta_vals::Vector{Float64})
    cevs = Float64[]
    iqr_d_income = Float64[]
    iqr_d_wealth = Float64[]
    iqr_d_c = Float64[]
    labels = String[]

    pe_baseline = copy(pe_base)
    pe_baseline[16] = 0.0  # θ = 0 baseline
    v_base = mean(valfun(pe_baseline).v)

    for thet in theta_vals
        pe_test = copy(pe_base)
        pe_test[16] = thet
        println("\n--- θ = $(thet) ---")
        result = valfun(pe_test)

        # === 1. Welfare: CEV ===
        v_base = vec(mean(valfun(pe_baseline).v, dims=(1,2,3)))
        v_post = vec(mean(result.v, dims=(1,2,3)))
        cev = compute_cev(v_base, v_post, pe_test)
        push!(cevs, cev)

        # === 2. Simulation + Moments (get dispersion) ===
        simdata = simmodel(result)
        moms, _, _, _, iqr_income, iqr_wealth, iqr_c = makemoments(simdata, pe_test)

        push!(iqr_d_income, iqr_income)
        push!(iqr_d_wealth, iqr_wealth)
        push!(iqr_d_c, iqr_c)
        
        push!(labels, "θ = $(round(thet, digits=2))")

        # === 3. Ergodic Distribution ===
        dist = compute_ergodic(result)
        plot_ergodic_a(dist, result.g)
        plot_ergodic_d(dist, result.g)
        plot_joint_heatmap(dist, result.g)
    end

    return labels, cevs, iqr_d_income, iqr_d_wealth, iqr_d_c
end


function run_batch()
    pe = ptrue(sz.nop)
    thetas = 0.0:0.5:1.0
    labels, cevs, iqr_income, iqr_wealth, iqr_c = batch_welfare_and_dispersion(pe, collect(thetas))

    println("\n--- Welfare Comparison (CEV %) ---")
    for (l, c) in zip(labels, cevs)
        println("$l ⇒ $(round(c, digits=2))%")
    end

    println("\n--- IQR: Durable Ratios ---")
    for i in 1:length(labels)
        println("$(labels[i]) → IQR[d/y] = $(round(iqr_income[i], digits=2)), IQR[d/w] = $(round(iqr_wealth[i], digits=2)), IQR[d/c] = $(round(iqr_c[i], digits=2))")
    end
    
end
