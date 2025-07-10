function utility_noadjust(grids::NamedTuple, pea::Vector{Float64})
    a, d, ap, dp, e, y = grids.a, grids.d, grids.ap, grids.dp, grids.ex, grids.y
    beta, delta, nu, gamma, w, chi, pd, tau, h = pea[1], pea[2], pea[5], pea[6], pea[8], pea[9], pea[10], pea[12], pea[13]
    rr = (1 / beta) - 1

    d_next = (1 - delta * (1 - chi)) .* d
    iid_map = [argmin(abs.(dp .- d_next[id])) for id in 1:sz.nd]

    util = zeros(sz.ne, sz.ny, sz.na, sz.nd, sz.npa, sz.npd)
    penalty_count = 0

    Threads.@threads for iia in 1:sz.npa
        Threads.@threads for id in 1:sz.nd
            Threads.@threads for ia in 1:sz.na
                Threads.@threads for iy in 1:sz.ny
                    Threads.@threads for ie in 1:sz.ne
                        d_next_val = d_next[id]
                        c = y[iy] * w * h * (1 - tau) + e[ie] * a[ia] * (1 + rr) - e[ie] * ap[iia]

                        if c > 0 && d_next_val > 0
                            util[ie, iy, ia, id, iia, iid_map[id]] = ((c^nu * d_next_val^(1 - nu))^(1 - gamma)) / (1 - gamma)
                        else
                            util[ie, iy, ia, id, iia, iid_map[id]] = -1e10
                            penalty_count += 1
                        end
                    end
                end
            end
        end
    end

    println("Number of penalized states NA: ", penalty_count)
    println("Share of penalized states NA: ", penalty_count / length(util))
    return util
end
