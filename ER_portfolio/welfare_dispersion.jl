function batch_welfare_and_dispersion(pe_base::Vector{Float64}, theta_vals::Vector{Float64})
    # Prepare output containers
    cevs = Float64[]
    iqr_d_wealth = Float64[]
    iqr_d_c = Float64[]
    adj_rates = Float64[]
    avg_durables = Float64[]    
    iqr_d = Float64[]
    p90_10_d_wealth = Float64[]
    p90_10_d_c = Float64[]
    p90_10_d = Float64[]
    mu_d_wealth = Float64[]
    mu_d_c = Float64[]
    mu_a = Float64[]  # Mean assets
    mu_gap = Float64[]  # Mean durables
    labels = String[]

    cev_map = Dict{Float64, Tuple{Array{Float64,4}, Array{Float64,4}}}()

    # === θ = 0 Baseline ===
    θ0 = 0.0
    pe_baseline = copy(pe_base)
    pe_baseline[16] = θ0
    result0 = valfun(pe_baseline)
    v_base = vec(mean(result0.v, dims=(1,2,3)))

    # Compute CEV relative to itself (0)
    push!(cevs, 0.0)  # baseline

    simdata0 = simmodel(result0)
    moms0, _, _, _, IQR_d_val, p90_10_d_val = makemoments(simdata0, pe_baseline)

    push!(iqr_d_wealth, moms0[13])
    push!(iqr_d_c, moms0[14])
    push!(iqr_d, IQR_d_val)
    push!(p90_10_d_wealth, moms0[15])
    push!(p90_10_d_c, moms0[16])
    push!(p90_10_d, p90_10_d_val)
    push!(avg_durables, moms0[1])
    push!(adj_rates, moms0[12])
    push!(mu_d_wealth, moms0[7])
    push!(mu_d_c, moms0[8])
    push!(mu_a, moms0[3])  # Mean assets
    push!(mu_gap, moms0[9])  # Mean durables
    push!(labels, "θ = 0.0")

    dist0 = compute_ergodic(result0)
 

    θ_folder0 = "Output/Ergodic/theta0.0"
    mkpath(θ_folder0)
    plot_ergodic_a(dist0, result0.g, θ_folder0)
    plot_ergodic_d(dist0, result0.g, θ_folder0)
    plot_joint_heatmap(dist0, result0.g, θ_folder0)
    plot_total_marginals(dist0, result0.g, θ_folder0)
    plot_adjustment_by_income(result0.adjustment_indicator, dist0, θ_folder0)

    # === Loop over remaining θs (excluding 0) ===
    for thet in theta_vals
        if isapprox(thet, 0.0); continue; end  # skip θ = 0 again

        pe_test = copy(pe_base)
        pe_test[16] = thet
        θ_label = "theta$(round(thet, digits=2))"
        θ_folder = "Output/Ergodic/$(θ_label)"
        mkpath(θ_folder)
        println("\n--- θ = $(thet) ---")
        result = valfun(pe_test)

        v_post = vec(mean(result.v, dims=(1,2,3)))
        cev = compute_cev(v_base, v_post, pe_test)
        push!(cevs, cev)

        simdata = simmodel(result)
        moms, _, _, _, IQR_d_val, p90_10_d_val = makemoments(simdata, pe_test)

        push!(iqr_d_wealth, moms[13])
        push!(iqr_d_c, moms[14])
        push!(iqr_d, IQR_d_val)
        push!(p90_10_d_wealth, moms[15])
        push!(p90_10_d_c, moms[16])
        push!(p90_10_d, p90_10_d_val)
        push!(avg_durables, moms[1])
        push!(adj_rates, moms[12])
        push!(mu_d_wealth, moms[7])
        push!(mu_d_c, moms[8])
        push!(mu_a, moms[3])  # Mean assets
        push!(mu_gap, moms[9])  # Mean durables
        push!(labels, "θ = $(round(thet, digits=2))")

        dist = compute_ergodic(result)
        cev_array = compute_cev_array(result.v, v_base, pe_test)
        cev_map[thet] = (cev_array, dist)

        plot_ergodic_a(dist, result.g, θ_folder)
        plot_ergodic_d(dist, result.g, θ_folder)
        plot_joint_heatmap(dist, result.g, θ_folder)
        plot_total_marginals(dist, result.g, θ_folder)
        plot_adjustment_by_income(result.adjustment_indicator, dist, θ_folder)
    end

    return (
        labels, cevs, iqr_d_wealth, iqr_d_c, adj_rates, avg_durables, cev_map,
        iqr_d, p90_10_d_wealth, p90_10_d_c, p90_10_d, mu_d_wealth, mu_d_c, mu_a, mu_gap
    )
end




"""
    plot_cev_distributions_all(cev_map)

Plots overlaid CEV welfare distributions for all θ cases.
"""
function plot_cev_distributions_all(cev_map::Dict{Float64, Tuple{Array{Float64,4}, Array{Float64,4}}})
    plt = plot(title="Welfare Distribution (CEV) by Currency Regime", 
               xlabel="CEV (%)", ylabel="Density", legend=:topright)

    # Sort by θ value ascending
    for (θ, (cev_array, dist)) in sort(collect(cev_map); by=x -> x[1])
        cev_vec = cev_array[:]
        weights = dist[:]
        weights ./= sum(weights)  # normalize weights

        histogram!(plt, cev_vec .* 100,  # convert to %
                   weights=weights,
                   bins=60, normalize=true,
                   alpha=0.5, lw=2,
                   label="θ = $(round(θ, digits=2))")
    end

    savefig(plt, "Output/Ergodic/cev_distribution_all.png")
end


"""
Saves a LaTeX-formatted table of welfare and dispersion statistics, with θ as columns.
"""
function save_latex_welfare_table(
    labels, cevs, 
    iqr_d_wealth, iqr_d_c, iqr_d, 
    p90_10_d_wealth, p90_10_d_c, p90_10_d, 
    mu_d_wealth, mu_d_c, 
    adj_rates, avg_durables, mu_a, mu_gap,
    filepath::String
)
    open(filepath, "w") do io
        # Extract just the numeric θs for columns
        thetas = [split(lbl, "=")[2] |> x -> strip(x) for lbl in labels]

        println(io, "\\begin{tabular}{l" * "c"^length(thetas) * "}")
        println(io, "\\toprule")
        println(io, "& " * join(["\$\\theta = $θ\$" for θ in thetas], " & ") * " \\\\")
        println(io, "\\midrule")

        function print_row(label::String, data::Vector{Float64}; percentage=false)
            row = percentage ? [round(x * 100, digits=2) for x in data] : [round(x, digits=2) for x in data]
            println(io, label * " & " * join(row, " & ") * " \\\\")
        end

        print_row("CEV (\\%)", cevs)
        print_row("IQR[d/w]", iqr_d_wealth)
        print_row("IQR[d/c]", iqr_d_c)
        print_row("IQR[d]", iqr_d)
        print_row("p90-10[d/w]", p90_10_d_wealth)
        print_row("p90-10[d/c]", p90_10_d_c)
        print_row("p90-10[d]", p90_10_d)
        print_row("Avg.[d/w]", mu_d_wealth)
        print_row("Avg. [d, c]", mu_d_c)
        print_row("Adj. Rate (\\%)", adj_rates; percentage=true)
        print_row("Change Durables", avg_durables)
        print_row("Change Assets", mu_a)
        print_row("Mean Gap", mu_gap)

        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
end


"""
    run_batch()

Example runner to execute batch welfare and dispersion analysis,
plot the CEV distributions, and save a LaTeX table.
"""
function run_batch()
    pe = ptrue(sz.nop)
    thetas = 0.25:0.25:1.0

    labels, cevs, iqr_d_wealth, iqr_d_c, adj_rates, avg_durables, cev_map,
        iqr_d, p90_10_d_wealth, p90_10_d_c, p90_10_d, mu_d_wealth, mu_d_c, mu_a, mu_gap = batch_welfare_and_dispersion(pe, collect(thetas))

    println("\n--- Welfare Gains (vs θ = 0 Baseline) ---")
    for (l, c) in zip(labels, cevs)
        println("$l ⇒ $(round(c, digits=2))%")
    end

    # Plot all CEV distributions together
    plot_cev_distributions_all(cev_map)

    # Plot each CEV distribution separately
    for θ in sort(collect(keys(cev_map)))
        cev_array, dist = cev_map[θ]
        plot_cev_distribution(cev_array, dist,θ)
        savefig("Output/Ergodic/cev_distribution_theta$(round(θ,digits=2)).png")
    end

    # Export LaTeX table
    save_latex_welfare_table(labels, cevs, 
    iqr_d_wealth, iqr_d_c, iqr_d, 
    p90_10_d_wealth, p90_10_d_c, p90_10_d, 
    mu_d_wealth, mu_d_c, 
    adj_rates, avg_durables, mu_a, mu_gap, "Output/Ergodic/welfare_table.tex")
end
