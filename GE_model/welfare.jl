function compute_cev(v_pre, v_post, ppp::Vector{Float64})
    nu    = ppp[5]
    gamma = ppp[6]

    mean_v_pre  = isa(v_pre, AbstractVector) ? mean(v_pre) : v_pre
    mean_v_post = isa(v_post, AbstractVector) ? mean(v_post) : v_post

    ratio = mean_v_post / mean_v_pre
    cev   = (ratio)^(1 / ((1 - gamma) * (1 - nu))) - 1
    return cev * 100
end

function welfare_comparison(pe_base::Vector{Float64}, theta_vals::Vector{Float64})
    cevs = Float64[]
    labels = String[]

    # Baseline (thet = 0)
    pe_baseline = copy(pe_base)
    pe_baseline[16] = 0.0
    ans_base = valfun(pe_baseline)
    v_base = mean(ans_base.v)

    for thet in theta_vals
        pe_test = copy(pe_base)
        pe_test[16] = thet

        ans_test = valfun(pe_test)
        v_test = mean(ans_test.v)

        cev = compute_cev(v_base, v_test, pe_test)
        push!(cevs, cev)
        push!(labels, "θ = $(round(thet, digits=2))")
    end

    return labels, cevs
end


function run_welfare_analysis()
    pea = ptrue(sz.nop)
    theta_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
    labels, cevs = welfare_comparison(pea, theta_vals)

    for (lbl, cev) in zip(labels, cevs)
        println("$lbl => Welfare gain: $(round(cev, digits=3))%")
    end
end
