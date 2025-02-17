function maxbellman_noadjust(queuelong::Array{Float64}, util::Array{Float64})
    vnew = zeros(sz.np, sz.ne, sz.na, sz.nd)
    gidx = dtp.Ipol(Int.(zeros(sz.np, sz.ne, sz.na, sz.nd)), Int.(zeros(sz.np, sz.ne, sz.na, sz.nd)))
    beta = pea[1]

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    vstar = -Inf
                    astar = 0
                    for iia in 1:sz.npa
                        bellman = util[ip, ie, ia, id, iia] + beta * queuelong[ip, ie, iia, id]
                        if bellman > vstar
                            vstar = bellman
                            astar = iia
                        end
                    end
                    vnew[ip, ie, ia, id] = vstar
                    gidx.a[ip, ie, ia, id] = astar
                end
            end
        end
    end

    # Update gidx.d outside the loops
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    gidx.d[ip, ie, ia, id] = id
                end
            end
        end
    end

    return vnew, gidx
end