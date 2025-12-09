"""
    welfare_disaster(pe; pi_annual, kappa_e)

Baseline vs. disaster counterfactual with SAME parameter vector `pe` but
different exchange-rate process.

Uses the same decomposition as `welfare_full_summary`, but:
- A: baseline grids (no disaster)
- B: disaster grids with (pi_annual, kappa_e)
"""
function welfare_disaster(pe::Vector{Float64};
                          pi_annual::Float64,
                          kappa_e::Float64)

    nu    = pe[5]
    gamma = pe[6]

    # We'll temporarily switch grid_builder
    global grid_builder

    # --- Baseline A ---
    grid_builder = p -> makegrids(p)
    ansA = valfun(pe)
    μA5  = compute_ergodic(ansA)                # [ie,iy,id,iaa,ia]
    μA   = permutedims(μA5, (1,2,4,5,3))        # [ie,iy,iaa,ia,id]

    # --- Disaster B ---
    grid_builder = p -> makegrids_disaster(p;
                                           pi_annual = pi_annual,
                                           kappa_e   = kappa_e)
    ansB = valfun(pe)
    μB5  = compute_ergodic(ansB)
    μB   = permutedims(μB5, (1,2,4,5,3))

    # Restore baseline builder for safety
    grid_builder = p -> makegrids(p)

    # Welfare objects
    wAA = sum(ansA.v .* μA)   # A under A grids
    wBB = sum(ansB.v .* μB)   # B under B grids
    wBA = sum(ansB.v .* μA)   # B value on A's distribution
    wAB = sum(ansA.v .* μB)   # A value on B's distribution

    # CEVs (composite)
    cev_BA = (wBB / wAA)^(1/(1 - gamma)) - 1
    cev_AB = (wAA / wBB)^(1/(1 - gamma)) - 1

    # Welfare changes (percent)
    keepDistAB = (wBA - wAA) / abs(wAA) * 100
    keepDistBA = (wAB - wBB) / abs(wBB) * 100
    acrossSS   = (wBB - wAA) / abs(wAA) * 100

    λ_composite = (wAA / wBB)^(1/(1 - gamma)) - 1
    λ_c_only    = (wAA / wBB)^(1/(nu * (1 - gamma))) - 1

    return (wAA=wAA, wBA=wBA, wBB=wBB, wAB=wAB,
            cev_BA=cev_BA, cev_AB=cev_AB,
            λ_composite=λ_composite, λ_c_only=λ_c_only,
            keepDistAB=keepDistAB, keepDistBA=keepDistBA, acrossSS=acrossSS)
end
