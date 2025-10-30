function utility_noadjust(grids::NamedTuple, pea::Vector{Float64})
    a, aa, d, ap, aap, dp, e, y =
        grids.a, grids.aa, grids.d, grids.ap, grids.aap, grids.dp, grids.ex, grids.y

    beta, delta, nu, gamma = pea[1], pea[2], pea[5], pea[6]
    w, pd, kappa, tau, h      = pea[8], pea[10], pea[11], pea[12], pea[13]
    rr = (1 / beta) - 1
    rr_star = pea[9]

    d_next = (1 - delta ) .* d
    iid_map = [argmin(abs.(dp .- d_next[id])) for id in 1:sz.nd]

    # util[ie,iy,iaa,ia,id, iiaa,iia, iid_map[id]]
    util = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd, sz.npa, sz.npa, sz.npd)

    Threads.@threads for iia in 1:sz.npa
         for iiaa in 1:sz.npa
             for id in 1:sz.nd
                 for ia in 1:sz.na
                     for iaa in 1:sz.na
                        for iy in 1:sz.ny
                            for ie in 1:sz.ne
                                E  = e[ie]; Y = y[iy]
                                aa_now = aa[iaa]
                                a_now  = a[ia]

                                carry = (1 + rr) * aa_now + (1+rr_star)*(E*a_now)

                                aa_next = aap[iiaa]
                                a_next  = ap[iia]
                                next_pay = aa_next + E*a_next
                                dollar_cost = kappa*(E * a_next)

                                inc = Y*w*h*(1 - tau) + carry
                                c = inc - next_pay  - dollar_cost

                                idd = iid_map[id]
                                dkeep = d_next[id]

                                if c > 0 && dkeep > 0 
                                    util[ie, iy, iaa, ia, id, iiaa, iia, idd] =
                                        ((c^nu * dkeep^(1 - nu))^(1 - gamma)) / (1 - gamma)
                                else
                                    util[ie, iy, iaa, ia, id, iiaa, iia, idd] = -1e10
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
