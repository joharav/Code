"""
    welfare_disaster(pe; pi_annual, kappa_e_log, periods_per_year=4)

A = baseline grids
B = disaster grids (same pe, different exogenous process)
"""
function welfare_disaster(pe::Vector{Float64};
                          pi_annual::Float64,
                          kappa_e_log::Float64,
                          periods_per_year::Int = 4)

    gridA = makegrids

    gridB = function (q::Vector{Float64})
        g_dis, _meta = makegrids_disaster(q;
            pi_annual=pi_annual,
            kappa_e_log=kappa_e_log,
            periods_per_year=periods_per_year
        )
        return g_dis
    end

    res = welfare_summary(pe, pe; gridA=gridA, gridB=gridB)
    return res
end
