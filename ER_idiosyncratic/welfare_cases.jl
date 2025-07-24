using Roots, PrettyTables


"""
    solve_lambda_transfer(pe_other, μ_base, w_base, ppp)

Solves for λ that makes welfare in the `pe_other` policy regime
equal to `w_base`, using `μ_base` as the distribution.
"""
function solve_lambda_transfer(pe_other, μ_base, w_base, ppp)
    function welfare_in_other_with_transfer(λ)
        ans_new = valfun(pe_other; λ=λ)
        return sum(ans_new.v .* μ_base)
    end

    λ_star = find_zero(λ -> welfare_in_other_with_transfer(λ) - w_base,
                       (0.0, 1.0), Bisection())
    return λ_star
end



"""
    welfare_full_summary(pe_A, pe_B, ppp)

Computes:
- Baseline welfare in A and B
- Welfare from switching regimes (fixed distributions)
- Transfers λ needed to compensate
- PrettyTables report and NamedTuple summary
"""
function welfare_full_summary(pe_A::Vector{Float64}, pe_B::Vector{Float64}, ppp::Vector{Float64})
    γ = ppp[6]

    # Solve baseline (λ = 0)
    ansA = valfun(pe_A; λ=0.0)
    ansB = valfun(pe_B; λ=0.0)
    μA = compute_ergodic(ansA)
    μB = compute_ergodic(ansB)

    # Welfare in each case
    wAA = sum(ansA.v .* μA)    # (1) Welfare in A under A
    wBB = sum(ansB.v .* μB)    # (4) Welfare in B under B
    wBA = sum(ansB.v .* μA)    # (2) Welfare in B under A's distribution
    wAB = sum(ansA.v .* μB)    # (5) Welfare in A under B's distribution

    # CEV comparisons (aggregate)
    cev_BA = compute_cev(vec(ansA.v), vec(ansB.v), ppp)
    cev_AB = compute_cev(vec(ansB.v), vec(ansA.v), ppp)

    # Solve for compensating transfers (λ)
    λ_BA = solve_lambda_transfer(pe_B, μA, wAA, ppp)  # (3)
    λ_AB = solve_lambda_transfer(pe_A, μB, wBB, ppp)  # (6)

    # Solve model again under λ transfers for record
    ansBλ = valfun(pe_B; λ=λ_BA)
    ansAλ = valfun(pe_A; λ=λ_AB)
    wBλA = sum(ansBλ.v .* μA)  # should ≈ wAA
    wAλB = sum(ansAλ.v .* μB)  # should ≈ wBB

    # Pretty Print
    using PrettyTables
    println("\n📊 Full Welfare Summary (γ = $γ):\n")
    data = [
        ("(1) Baseline A",         wAA, "", ""),
        ("(2) Counterfactual B in A", wBA, "CEV = $(round(cev_BA, digits=2))%", ""),
        ("(3) Transfer in A to restore welfare", wBλA, "", "λ = $(round(λ_BA * 100, digits=2))%"),
        ("(4) Baseline B",         wBB, "", ""),
        ("(5) Counterfactual A in B", wAB, "CEV = $(round(cev_AB, digits=2))%", ""),
        ("(6) Transfer in B to restore welfare", wAλB, "", "λ = $(round(λ_AB * 100, digits=2))%"),
    ]
    pretty_table(data, header=["Case", "Welfare", "CEV", "Transfer λ"])

    return (
        wAA = wAA, wBA = wBA, λ_BA = λ_BA, cev_BA = cev_BA,
        wBB = wBB, wAB = wAB, λ_AB = λ_AB, cev_AB = cev_AB,
        wBλA = wBλA, wAλB = wAλB,
        μA = μA, μB = μB,
        ansA = ansA, ansB = ansB,
        ansBλ = ansBλ, ansAλ = ansAλ
    )
end


results = welfare_full_summary(pe_A, pe_B, ppp)

println("\n→ Case 2: CEV loss in A if switch to B: ", round(results.cev_BA, digits=2), "%")
println("→ Case 3: λ transfer needed in B to restore A's welfare: ", round(results.λ_BA * 100, digits=2), "%")

println("\n→ Case 5: CEV gain in B if switch to A: ", round(results.cev_AB, digits=2), "%")
println("→ Case 6: λ transfer needed in A to restore B's welfare: ", round(results.λ_AB * 100, digits=2), "%")
