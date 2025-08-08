function makemoments_target(simdata::NamedTuple, pea::Vector{Float64})
    beta  = pea[1]
    w     = pea[8]
    pd    = pea[10]
    theta = pea[16]
    rr    = 1/beta - 1

    # Extract variables (no adjustment indicators needed)
    a = simdata.a[sz.burnin:end, :]
    d = simdata.d[sz.burnin:end, :]
    y = simdata.y[sz.burnin:end, :]
    ex = simdata.ex[sz.burnin:end, :]

    # Effective wealth
    a_effective = theta * ex .* a + (1 - theta) .* a

    # Moments
    mu_d  = mean(vec(d))
    var_d = var(vec(d))
    mu_a  = mean(vec(a))
    var_a = var(vec(a))

    ratio_d_income = vec(pd .* ex .* d) ./ vec(w .* y .+ a_effective .* (1 + rr))
    ratio_d_wealth = vec(pd .* ex .* d) ./ vec(a_effective .* (1 + rr) .+ pd .* ex .* d)

    mu_d_income = mean(ratio_d_income)
    mu_d_wealth = mean(ratio_d_wealth)

    # Dispersion
    disp_d_income = compute_dispersion(ratio_d_income)
    disp_d_wealth = compute_dispersion(ratio_d_wealth)
    IQR_d_income  = disp_d_income[2]
    IQR_d_wealth  = disp_d_wealth[2]
    p90_10_d_income = disp_d_income[3]
    p90_10_d_wealth = disp_d_wealth[3]

    # Assemble target moments
    outmoms = [
        mu_d, var_d, mu_a, var_a,
        mu_d_income, mu_d_wealth,
        IQR_d_income, IQR_d_wealth,
        p90_10_d_income, p90_10_d_wealth
    ]

    return outmoms::Vector{Float64}
end
