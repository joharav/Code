function utility_noadjust(grids::NamedTuple, pea::Vector{Float64})
    a, aa, d, ap, aap, dp, e, y =
        grids.a, grids.aa, grids.d, grids.ap, grids.aap, grids.dp, grids.ex, grids.y

    beta, delta, nu, gamma = pea[1], pea[2], pea[5], pea[6]
    w, chi, pd, tau, h      = pea[8], pea[9], pea[10], pea[12], pea[13]
    rr = (1 / beta) - 1

    d_next = (1 - delta * (1 - chi)) .* d
    iid_map = [argmin(abs.(dp .- d_next[id])) for id in 1:sz.nd]

    # util[ie,iy,iaa,ia,id, iiaa,iia, iid_map[id]]
    util = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd, sz.npa, sz.npa, sz.npd)
    penalty_count = 0

    Threads.@threads for iia in 1:sz.npa
        Threads.@threads for iiaa in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                Threads.@threads for ia in 1:sz.na
                    Threads.@threads for iaa in 1:sz.na
                        Threads.@threads for iy in 1:sz.ny
                            Threads.@threads for ie in 1:sz.ne
                                E  = e[ie]; Y = y[iy]
                                aa_now = aa[iaa]
                                a_now  = a[ia]

                                carry = (1 + rr) * (aa_now + E*a_now)

                                aa_next = aap[iiaa]
                                a_next  = ap[iia]
                                next_pay = aa_next + E*a_next

                                maintenance = E * pd * delta * chi * d[id]
                                inc = Y*w*h*(1 - tau) + carry
                                c = inc - next_pay - maintenance

                                idd = iid_map[id]
                                dkeep = d_next[id]

                                if c > 0 && dkeep > 0 
                                    util[ie, iy, iaa, ia, id, iiaa, iia, idd] =
                                        ((c^nu * dkeep^(1 - nu))^(1 - gamma)) / (1 - gamma)
                                else
                                    util[ie, iy, iaa, ia, id, iiaa, iia, idd] = -1e10
                                    penalty_count += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if settings.verbose
        println("Number of penalized states NA: ", penalty_count)
        println("Share of penalized states NA: ", penalty_count / length(util))
    end
    return util
end
