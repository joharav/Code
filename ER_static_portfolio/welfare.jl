using Statistics

# Utility mapping: CRRA over composite consumption c with u(c)=c^(1-γ)/(1-γ)
# If your value function is already lifetime utility in those units, CEV is:
# λ = (W_B/W_A)^(1/(1-γ)) - 1
@inline function cev_from_W(WA::Float64, WB::Float64, γ::Float64)
    if γ == 1.0
        return exp(WB - WA) - 1.0
    else
        return (WB / WA)^(1.0 / (1.0 - γ)) - 1.0
    end
end

"""
    welfare_summary(peA, peB; gridA=makegrids, gridB=makegrids)

Solve two economies, compute ergodic distributions, and return
(wAA, wBA, wBB, wAB) plus CEVs and decompositions.

All objects are 4D (post-portfolio dynamic states).
"""
function welfare_summary(peA::Vector{Float64}, peB::Vector{Float64};
                         gridA::Function = makegrids,
                         gridB::Function = makegrids)

    γ = peA[6]  # keep consistent with your parameter vector

    ansA = valfun(peA)
    ansB = valfun(peB)

    μA = compute_ergodic(ansA)   # same shape as ansA.v
    μB = compute_ergodic(ansB)

    @assert size(μA) == size(ansA.v)
    @assert size(μB) == size(ansB.v)

    wAA = sum(ansA.v .* μA)
    wBA = sum(ansB.v .* μA)
    wBB = sum(ansB.v .* μB)
    wAB = sum(ansA.v .* μB)

    cev_BA = cev_from_W(wAA, wBB, γ)        # switch A→B in steady state
    cev_AB = cev_from_W(wBB, wAA, γ)

    acrossSS   = (wBB - wAA) / abs(wAA) * 100

    # “keep distribution” effects as CEV objects (not percent change in W)
    cev_keepDist_AB = cev_from_W(wAA, wBA, γ)  # evaluate B under μA
    cev_keepDist_BA = cev_from_W(wBB, wAB, γ)  # evaluate A under μB

    return (wAA=wAA, wBA=wBA, wBB=wBB, wAB=wAB,
            cev_BA=cev_BA, cev_AB=cev_AB,
            cev_keepDist_AB=cev_keepDist_AB,
            cev_keepDist_BA=cev_keepDist_BA, acrossSS=acrossSS)
end
