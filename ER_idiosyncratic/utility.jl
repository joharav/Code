function utility(grids::NamedTuple, pea::Vector{Float64}; λ::Float64 = 0.0)
    a, d, ap, dp, e, y = grids.a, grids.d, grids.ap, grids.dp, grids.ex, grids.y
  
    beta, delta, nu, gamma = pea[1], pea[2], pea[5], pea[6]
    f, w, pd, ft, tau, h    = pea[7], pea[8], pea[10], pea[11], pea[12], pea[13]
    rr = (1 / beta) - 1
    theta = pea[16]
    util = zeros(sz.ne, sz.ny, sz.na, sz.nd, sz.npa, sz.npd)
    penalty_count = 0

    Threads.@threads for iid in 1:sz.npd
        Threads.@threads for iia in 1:sz.npa
            Threads.@threads for id in 1:sz.nd
                Threads.@threads for ia in 1:sz.na
                    Threads.@threads for iy in 1:sz.ny
                        Threads.@threads for ie in 1:sz.ne
                            a_effective = theta * e[ie] * a[ia] + (1 - theta) * a[ia]
                            a_eff_prime = theta * e[ie] * ap[iia] + (1 - theta) * ap[iia]
                            income = y[iy] * w * h * (1 - tau) + a_effective * (1 + rr)
                            sale = e[ie] * pd * (1 - f) * (1 - delta) * d[id]
                            cost = e[ie] * pd * dp[iid]
                            timecost = y[iy] * w * h * ft
                            c = income + sale - cost - a_eff_prime - timecost

                            if c > 0 && dp[iid] > 0 && (1 + λ) > 0
                                c_eff = c * (1 + λ)
                                util[ie, iy, ia, id, iia, iid] = ((c_eff^nu * dp[iid]^(1 - nu))^(1 - gamma)) / (1 - gamma)
                            else
                                util[ie, iy, ia, id, iia, iid] = -1e10
                                penalty_count += 1
                            end
                        end
                    end
                end
            end
        end
    end

    penalty_share = penalty_count / (sz.ne * sz.ny * sz.na * sz.nd * sz.npa * sz.npd)
    if settings.verbose
        println("Number of penalized states ADJUST: ", penalty_count)
        println("Share of penalized states ADJUST: ", penalty_share)
    end



    return util
end
