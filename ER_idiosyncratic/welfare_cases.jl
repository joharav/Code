"""
    welfare_full_summary(pe_A, pe_B)

Computes:
- Baseline welfare in A and B
- Welfare from switching regimes (fixed distributions)
- Transfers λ needed to compensate
- PrettyTables report and NamedTuple summary
"""
function welfare_full_summary(pe_A::Vector{Float64}, pe_B::Vector{Float64})
    nu = pe_A[5]
    gamma = pe_A[6]
    # Solve baseline 
    ansA = valfun(pe_A)
    ansB = valfun(pe_B)
    μA = compute_ergodic(ansA)
    μB = compute_ergodic(ansB)
    simB=simmodel(ansB)
    moms, _, _, _, _, _ = makemoments(simB, pe_B; shock=false)
    vol_c_with = moms[5]
    vol_x_with = moms[6]

    # #No dollars
    # pe_B_nod = copy(pe_B); pe_B_nod[16] = 0.0
    # ansB_nod = valfun(pe_B_nod)
    # simB_nod = simmodel(ansB_nod)
    # momsB_nod, _, _, _, _, _ = makemoments(simB_nod, pe_B_nod; shock=false)
    # vol_c_no = momsB_nod[5]
    # vol_x_no = momsB_nod[6]


    # Welfare in each case
    wAA = sum(ansA.v .* μA)    # (1) Welfare in A under A
    wBB = sum(ansB.v .* μB)    # (4) Welfare in B under B
    wBA = sum(ansB.v .* μA)    # (2) Welfare in B under A's distribution
    wAB = sum(ansA.v .* μB)    # (5) Welfare in A under B's distribution

    # CEV comparisons (aggregate)
  #  cev_BA = compute_cev(vec(ansA.v), vec(ansB.v), pe_A)
   # cev_AB = compute_cev(vec(ansB.v), vec(ansA.v), pe_A)

    # CEVs (composite-equivalent by default)
    cev_BA = (wBB / wAA)^(1/(1 - gamma)) - 1
    cev_AB = (wAA / wBB)^(1/(1 - gamma)) - 1


    # Welfare changes keeping distributions fixed
    keepDistAB = (wBA - wAA) / abs(wAA) * 100   # % change A → B, keep A's μ
    keepDistBA = (wAB - wBB) / abs(wBB) * 100   # % change B → A, keep B's μ
    acrossSS   = (wBB - wAA) / abs(wAA) * 100  # across steady states

    # Compensating transfer (choose one; label in slides!)
    λ_composite = (wAA / wBB)^(1/(1 - gamma)) - 1
    # or, if you truly want c-only:
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
    # println("Comp. transfer λ (c-only):    $(round(λ_c_only*100, digits=2))%")

    return (wAA=wAA, wBA=wBA, wBB=wBB, wAB=wAB,
    cev_BA=cev_BA, cev_AB=cev_AB,
    λ_composite=λ_composite, λ_c_only=λ_c_only,
    keepDistAB=keepDistAB, keepDistBA=keepDistBA, acrossSS=acrossSS, vol_c_with=vol_c_with, vol_x_with=vol_x_with) 
    
    #,
    #vol_c_with=vol_c_with, vol_x_with=vol_x_with,
    #vol_c_no=vol_c_no,     vol_x_no=vol_x_no
end
