# ============================================================
# welfare.jl  (FIXED)
# ============================================================

@inline function cev_from_W(WA::Float64, WB::Float64, γ::Float64)
    if γ == 1.0
        return exp(WB - WA) - 1.0
    else
        return (WB / WA)^(1.0 / (1.0 - γ)) - 1.0
    end
end

function welfare_summary(peA::Vector{Float64}, peB::Vector{Float64};
                         gridA::Function = makegrids,
                         gridB::Function = makegrids)

    γ = peA[6]

    # CRITICAL: pass the grid builders into valfun
    ansA = valfun(peA; grid_builder=gridA)
    ansB = valfun(peB; grid_builder=gridB)

    μA = compute_ergodic(ansA)
    μB = compute_ergodic(ansB)

    @assert size(μA) == size(ansA.v)
    @assert size(μB) == size(ansB.v)

    wAA = sum(ansA.v .* μA)
    wBA = sum(ansB.v .* μA)
    wBB = sum(ansB.v .* μB)
    wAB = sum(ansA.v .* μB)

    cev_BA = cev_from_W(wAA, wBB, γ)
    cev_AB = cev_from_W(wBB, wAA, γ)

    acrossSS = (wBB - wAA) / abs(wAA) * 100.0

    cev_keepDist_AB = cev_from_W(wAA, wBA, γ)
    cev_keepDist_BA = cev_from_W(wBB, wAB, γ)

    return (wAA=wAA, wBA=wBA, wBB=wBB, wAB=wAB,
            cev_BA=cev_BA, cev_AB=cev_AB,
            cev_keepDist_AB=cev_keepDist_AB,
            cev_keepDist_BA=cev_keepDist_BA,
            acrossSS=acrossSS)
end
