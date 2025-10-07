function utility(grids::NamedTuple, pea::Vector{Float64})
    a, aa, d, ap, aap, dp, e, y = grids.a, grids.aa, grids.d, grids.ap, grids.aap, grids.dp, grids.ex, grids.y
  
    beta, delta, nu, gamma = pea[1], pea[2], pea[5], pea[6]
    f, w, pd, kappa, tau, h    = pea[7], pea[8], pea[10], pea[11], pea[12], pea[13]
    rr = (1 / beta) - 1
    rr_star = pea[9]

    # util[ie,iy,iaa,ia,id, iiaa,iia,iid]
    util = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd, sz.npa, sz.npa, sz.npd)

    Threads.@threads for iid in 1:sz.npd
        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for iiaa in 1:sz.npa
                Threads.@threads for id in 1:sz.nd
                    Threads.@threads for ia in 1:sz.na
                        Threads.@threads for iaa in 1:sz.na
                            Threads.@threads for iy in 1:sz.ny
                                Threads.@threads for ie in 1:sz.ne
                                    # current holdings (local currency units)
                                    E  = e[ie]; Y = y[iy]
                                    aa_now = aa[iaa]
                                    a_now  = a[ia]

                                    # carryover resources (local currency): apply return to both assets
                                    carry = (1 + rr) * aa_now + (1+rr_star)*(E*a_now)

                                    # next holdings (local currency) chosen on policy grids
                                    aa_next = aap[iiaa]
                                    a_next  = ap[iia]
                                    next_pay = aa_next + E*a_next

                                    # income + durable adjust flows
                                    income = Y*w*h*(1 - tau) + carry
                                    sale_value = E * pd * (1 - f) * (1 - delta) * d[id]
                                    durable_purchase = E * pd * dp[iid]
                                    dollar_cost = kappa*(E * abs(a_next-a_now))

                                    c = income + sale_value - durable_purchase - next_pay - dollar_cost

                                    if c > 0 && dp[iid] > 0 
                                        util[ie, iy, iaa, ia, id, iiaa, iia, iid] =
                                            ((c^nu * dp[iid]^(1 - nu))^(1 - gamma)) / (1 - gamma)
                                    else
                                        util[ie, iy, iaa, ia, id, iiaa, iia, iid] = -1e10
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    if settings.verbose
        bad = count(u -> u == -1e10, util)
        println("Number penalized ... ", bad)
        println("Share penalized ... ", bad / length(util))
    end
    
    return util
end

