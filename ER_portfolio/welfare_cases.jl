"""
    welfare_full_summary(pe_A, pe_B)

Computes:
- Baseline welfare in A and B
- Welfare from switching regimes (fixed distributions)
- Transfers λ needed to compensate
- PrettyTables report and NamedTuple summary
"""
function welfare_full_summary(pe_A::Vector{Float64}, pe_B::Vector{Float64})
    nu    = pe_A[5]
    gamma = pe_A[6]

    # Solve baseline (λ = 0)
    ansA = valfun(pe_A)
    ansB = valfun(pe_B)

    μA5 = compute_ergodic(ansA)  # shape: [ie, iy, id, iaa, ia]
    μB5 = compute_ergodic(ansB)

    # Reorder to match v: [ie, iy, iaa, ia, id]
    μA = permutedims(μA5, (1,2,4,5,3))
    μB = permutedims(μB5, (1,2,4,5,3))

    # Welfare in each case
    wAA = sum(ansA.v .* μA)   # (1) Welfare in A under A
    wBB = sum(ansB.v .* μB)   # (4) Welfare in B under B
    wBA = sum(ansB.v .* μA)   # (2) Welfare in B under A's dist
    wAB = sum(ansA.v .* μB)   # (5) Welfare in A under B's dist

    # CEVs (composite-equivalent by default)
    cev_BA = (wBB / wAA)^(1/(1 - gamma)) - 1
    cev_AB = (wAA / wBB)^(1/(1 - gamma)) - 1

    # Welfare changes keeping distributions fixed
    keepDistAB = (wBA - wAA) / abs(wAA) * 100
    keepDistBA = (wAB - wBB) / abs(wBB) * 100
    acrossSS   = (wBB - wAA) / abs(wAA) * 100

    # Compensating transfer (composite)
    λ_composite = (wAA / wBB)^(1/(1 - gamma)) - 1
    # If you need c-only:
    λ_c_only    = (wAA / wBB)^(1/(nu * (1 - gamma))) - 1

    println("\n📊 Full Welfare Summary:\n")
    pretty_table([  "Baseline A (wAA)"               wAA    "";
                    "B on A's dist (wBA)"            wBA    "CEV (composite) = $(round(cev_BA*100, digits=2))%";
                    "Baseline B (wBB)"               wBB    "";
                    "A on B's dist (wAB)"            wAB    "CEV (composite) = $(round(cev_AB*100, digits=2))%"],
                 header=["Case","Welfare","Note"])

    println("\nAcross steady states (A→B): $(round(acrossSS, digits=2))%")
    println("A→B holding μ_A:            $(round(keepDistAB, digits=2))%")
    println("B→A holding μ_B:            $(round(keepDistBA, digits=2))%")
    println("Comp. transfer λ (composite): $(round(λ_composite*100, digits=2))%")

    return (wAA=wAA, wBA=wBA, wBB=wBB, wAB=wAB,
            cev_BA=cev_BA, cev_AB=cev_AB,
            λ_composite=λ_composite, λ_c_only=λ_c_only,
            keepDistAB=keepDistAB, keepDistBA=keepDistBA, acrossSS=acrossSS)
end
