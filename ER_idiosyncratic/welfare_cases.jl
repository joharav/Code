using Roots, PrettyTables


"""
    solve_lambda_transfer(pe_other, μ_base, w_base)

Solves for λ that makes welfare in the `pe_other` policy regime
equal to `w_base`, using `μ_base` as the distribution.
"""
function solve_lambda_transfer(pe_other, μ_base, w_base)
    function welfare_in_other_with_transfer(λ)
        ans_new = valfun(pe_other; λ=λ)
        return sum(ans_new.v .* μ_base)
    end

    λ_star = find_zero(λ -> welfare_in_other_with_transfer(λ) - w_base,
                       (-0.85, 1.0), Bisection())
    return λ_star
end



"""
    welfare_full_summary(pe_A, pe_B)

Computes:
- Baseline welfare in A and B
- Welfare from switching regimes (fixed distributions)
- Transfers λ needed to compensate
- PrettyTables report and NamedTuple summary
"""
function welfare_full_summary(pe_A::Vector{Float64}, pe_B::Vector{Float64})
    γ = pe_A[6]

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
    cev_BA = compute_cev(vec(ansA.v), vec(ansB.v), pe_A)
    cev_AB = compute_cev(vec(ansB.v), vec(ansA.v), pe_A)

   # Pretty Print
   println("\n📊 Full Welfare Summary (γ = $γ):\n")

   # --- FIX: Removed commas to create a 3-column Matrix instead of a 1-column Vector of Tuples ---
   data = [
       "(1) Baseline A"          wAA   "";
       "(2) Counterfactual B in A" wBA   "CEV = $(round(cev_BA, digits=2))%";
       "(3) Baseline B"          wBB   "";
       "(4) Counterfactual A in B" wAB   "CEV = $(round(cev_AB, digits=2))%"
   ]
   pretty_table(data, header=["Case", "Welfare", "CEV"])
# Correctly compares the change from baseline A to baseline B
println("\n→ Welfare Change from A to B: ", round((wBB - wAA) / abs(wAA) * 100, digits=2), "%")

# Correctly compares the change from baseline B to the counterfactual welfare in A
println("→ Welfare Change from A to B (keep ergodic distr): ", round((wBA - wAA) / abs(wAA) * 100, digits=2), "%")

    return (
        wAA = wAA, wBA = wBA, wBB = wBB, wAB = wAB,
        cev_BA = cev_BA, cev_AB = cev_AB
    )
end

pea = ptrue(sz.nop)
pe_A = pea  # Policy regime A
pe_B = copy(pe_A)
pe_B[16] = 0
results = welfare_full_summary(pe_A, pe_B)
