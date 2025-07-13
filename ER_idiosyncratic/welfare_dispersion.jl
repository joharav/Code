function batch_welfare_and_dispersion(pe_base::Vector{Float64}, theta_vals::Vector{Float64})
    cevs = Float64[]
    iqr_d_income = Float64[]
    iqr_d_wealth = Float64[]
    adj_rates = Float64[]
    mu_durables_change = Float64[]    
    iqr_d_c = Float64[]
    labels = String[]
    cev_map = Dict{Float64, Tuple{Array{Float64,4}, Array{Float64,4}}}()

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

        # === 2. Simulation + Moments ===
        simdata = simmodel(result)
        moms, _, _, _, iqr_income, iqr_wealth, iqr_c = makemoments(simdata, pe_test)

        push!(iqr_d_income, iqr_income)
        push!(iqr_d_wealth, iqr_wealth)
        push!(iqr_d_c, iqr_c)
        push!(adj_rates, moms[13])
        push!(mu_durables_change, moms[1])
        
        
        push!(labels, "θ = $(round(thet, digits=2))")

        # === 3. Ergodic Distribution ===
        dist = compute_ergodic(result)

        # === 4. CEV Distribution Array ===
        cev_array = compute_cev_array(result.v, v_base, pe_test)  # this must return 4D array of CEVs
        cev_map[thet] = (cev_array, dist)

        # Optional: Save ergodic plots per case
        plot_ergodic_a(dist, result.g)
        plot_ergodic_d(dist, result.g)
        plot_joint_heatmap(dist, result.g)
    end

    # Plot all CEV distributions together
    plot_cev_distributions_all(cev_map)

    return labels, cevs, iqr_d_income, iqr_d_wealth, iqr_d_c, adj_rates, avg_durables, cev_map
end



function plot_cev_distributions_all(cev_map::Dict{Float64, Tuple{Array{Float64,4}, Array{Float64,4}}})
    plt = plot(title="Welfare Distribution (CEV) by Currency Regime", 
               xlabel="CEV (%)", ylabel="Density", legend=:topright)

    for (θ, (cev_array, dist)) in sort(collect(cev_map); by=x -> x[1])
        cev_vec = cev_array[:]
        weights = dist[:]
        weights ./= sum(weights)

        histogram!(plt, cev_vec .* 100,  # CEV in percent
                   weights=weights,
                   bins=60, normalize=true,
                   alpha=0.5, lw=2,
                   label="θ = $(round(θ, digits=2))")
    end

    savefig(plt, "Output/Ergodic/cev_distribution_all.png")
end

function save_latex_welfare_table(labels, cevs, iqr_income, iqr_wealth, iqr_c, adj_rates, mu_durables_change, filepath::String)
    open(filepath, "w") do io
        println(io, "\\begin{tabular}{lcccccc}")
        println(io, "\\toprule")
        println(io, "\$\\theta\$ & CEV (\\%) & IQR[d/y] & IQR[d/w] & Adj. Rate & Avg. Δd & IQR[d/c] \\\\")
        println(io, "\\midrule")
        for i in 1:length(labels)
            θ = split(labels[i], "=")[2]
            cev = round(cevs[i], digits=2)
            iqry = round(iqr_income[i], digits=2)
            iqrw = round(iqr_wealth[i], digits=2)
            adj = round(adj_rates[i] * 100, digits=1)
            mud = round(mu_durables_change[i], digits=2)
            iqrc = round(iqr_c[i], digits=2)
            println(io, "\$\\theta = $θ\$ & $cev & $iqry & $iqrw & $adj\\% & $mud & $iqrc \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end



function run_batch()
    pe = ptrue(sz.nop)
    thetas = 0.0:0.5:1.0

    labels, cevs, iqr_income, iqr_wealth, iqr_c, adj_rates, mu_durables_change, cev_map = batch_welfare_and_dispersion(pe, collect(thetas))

    println("\n--- Welfare Gains (vs θ = 0 Baseline) ---")
    for (l, c) in zip(labels, cevs)
        println("$l ⇒ $(round(c, digits=2))%")
    end

    println("\n--- IQR: Durable Ratios ---")
    for i in 1:length(labels)
        println("$(labels[i]) → IQR[d/y] = $(round(iqr_income[i], digits=2)), IQR[d/w] = $(round(iqr_wealth[i], digits=2)), IQR[d/c] = $(round(iqr_c[i], digits=2))")
    end

    # === Plot all CEV distribution overlays ===
    plot_cev_distributions_all(cev_map)

    # === Plot each CEV distribution separately ===
    for θ in sort(collect(keys(cev_map)))
        cev_array, dist = cev_map[θ]
        plot_cev_distribution(cev_array, dist)
        savefig("Output/Ergodic/cev_distribution_theta$(θ).png")
    end

    # === Export LaTeX table ===
    save_latex_welfare_table(labels, cevs, iqr_income, iqr_wealth, iqr_c,adj_rates, avg_durables, "Output/Ergodic/welfare_table.tex")
end
