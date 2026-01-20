# ==========================================================================
# 4D MODEL: Combined VFI wrapper
# Solves both adjust and non-adjust regimes, then merges
# ==========================================================================

function valfun(pea::Vector{Float64})
    # Solve adjustment regime
    if settings.verbose
        println("\n=== Solving ADJUSTMENT regime ===")
    end
    adj_result = valfun_adjust(pea)
    
    # Solve non-adjustment regime
    if settings.verbose
        println("\n=== Solving NO-ADJUSTMENT regime ===")
    end
    nadj_result = valfun_noadjust(pea)
    
    # Merge: take max of values at each state
    # The merged value represents the optimal choice between adjusting or not
    v_merged = max.(adj_result.v, nadj_result.v)
    
    # For merged policy, choose from the winning regime at each state
    pol_merged = dtp.Pol(
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # w
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # d
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # s
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # c
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # aa
        zeros(sz.ne, sz.ny, sz.nw, sz.nd)    # a
    )
    
    for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
        if adj_result.v[ie, iy, iw, id] >= nadj_result.v[ie, iy, iw, id]
            pol_merged.w[ie, iy, iw, id] = adj_result.pol.w[ie, iy, iw, id]
            pol_merged.d[ie, iy, iw, id] = adj_result.pol.d[ie, iy, iw, id]
            pol_merged.s[ie, iy, iw, id] = adj_result.pol.s[ie, iy, iw, id]
            pol_merged.c[ie, iy, iw, id] = adj_result.pol.c[ie, iy, iw, id]
            pol_merged.aa[ie, iy, iw, id] = adj_result.pol.aa[ie, iy, iw, id]
            pol_merged.a[ie, iy, iw, id] = adj_result.pol.a[ie, iy, iw, id]
        else
            pol_merged.w[ie, iy, iw, id] = nadj_result.pol.w[ie, iy, iw, id]
            pol_merged.d[ie, iy, iw, id] = nadj_result.pol.d[ie, iy, iw, id]
            pol_merged.s[ie, iy, iw, id] = nadj_result.pol.s[ie, iy, iw, id]
            pol_merged.c[ie, iy, iw, id] = nadj_result.pol.c[ie, iy, iw, id]
            pol_merged.aa[ie, iy, iw, id] = nadj_result.pol.aa[ie, iy, iw, id]
            pol_merged.a[ie, iy, iw, id] = nadj_result.pol.a[ie, iy, iw, id]
        end
    end
    
    # Compute adjustment probability
    adj_prob = mean(adj_result.v .>= nadj_result.v)
    if settings.verbose
        println("\nAdjustment probability (value-based): ", round(adj_prob, digits=3))
    end
    
    return (
        v = v_merged,
        pol = pol_merged,
        g = adj_result.g,
        e = max(adj_result.e, nadj_result.e),
        adjust_result = adj_result,
        noadjust_result = nadj_result,
        pea = pea
    )
end
