function valfun(pea::Vector{Float64})
    noadjust_result = valfun_noadjust(pea)
    if settings.verbose
        println("Non-adjust: mean(v) = ", mean(noadjust_result.v))
        println("Non-adjust: min(v) = ", minimum(noadjust_result.v))
        println("Non-adjust: max(v) = ", maximum(noadjust_result.v))
    end

    adjust_result = valfun_adjust(pea)
    if settings.verbose
        println("Adjust: mean(v) = ", mean(adjust_result.v))
        println("Adjust: min(v) = ", minimum(adjust_result.v))
        println("Adjust: max(v) = ", maximum(adjust_result.v))
    end

    v = max.(adjust_result.v, noadjust_result.v)
    val_diff = adjust_result.v .- noadjust_result.v

    if settings.verbose
        println("Mean diff: ", mean(val_diff), " | Min: ", minimum(val_diff), " | Max: ", maximum(val_diff))
        println("Share where adjust gives lower value: ", sum(val_diff .< 0) / length(val_diff))
    end

    indicator_matrix = adjust_result.v .> noadjust_result.v
    adjust_flag = Float64.(indicator_matrix)

    if settings.verbose
        println("Pct adjusted policy chosen: ", 100*sum(indicator_matrix)/length(indicator_matrix))
    end

    # --- merge policies & indices (foreign a, local aa) ---
    gidx_a   = ifelse.(indicator_matrix, adjust_result.gidx.a,  noadjust_result.gidx.a)
    gidx_aa  = ifelse.(indicator_matrix, adjust_result.gidx.aa, noadjust_result.gidx.aa)
    gidx_d   = ifelse.(indicator_matrix, adjust_result.gidx.d,  noadjust_result.gidx.d)
    gidx = dtp.Ipol(gidx_a, gidx_aa, gidx_d)

    pol_a   = ifelse.(indicator_matrix, adjust_result.pol.a,  noadjust_result.pol.a)
    pol_aa  = ifelse.(indicator_matrix, adjust_result.pol.aa, noadjust_result.pol.aa)
    pol_d   = ifelse.(indicator_matrix, adjust_result.pol.d,  noadjust_result.pol.d)
    pol_c   = ifelse.(indicator_matrix, adjust_result.pol.c,  noadjust_result.pol.c)
    pol = dtp.Pol(pol_a, pol_aa, pol_d, pol_c)

    g = adjust_result.g
    e = adjust_result.e

    outtuple = (v=v, gidx=gidx, pol=pol, g=g, e=e,
                adjust_result=adjust_result, noadjust_result=noadjust_result,
                adjust_flag=adjust_flag, adjustment_indicator=indicator_matrix)

    if settings.verbose && adjust_result.e == 0 && noadjust_result.e == 0
        printstuff(outtuple)
        plotstuff(outtuple.v, gidx_a, gidx_d, outtuple.pol.c, g)
    end
    return outtuple
end
