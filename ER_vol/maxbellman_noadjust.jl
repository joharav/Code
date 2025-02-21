function maxbellman_noadjust(queuelong::Array{Float64}, util::Array{Float64})
    vnew = zeros(sz.ne, sz.na, sz.nd)
    gidx = dtp.Ipol(Int.(zeros(sz.ne, sz.na, sz.nd)), Int.(zeros(sz.ne, sz.na, sz.nd)))
    beta = pea[1]

    # Parallelize only on outermost loop
    Threads.@threads for id in 1:sz.nd
        for ia in 1:sz.na
            for ie in 1:sz.ne
                vstar = -1e20  # Avoid numerical instability
                astar = 0
                for iia in 1:sz.npa
                    bellman = util[ie, ia, id, iia] + beta * queuelong[ie, iia, id]
                    if bellman > vstar
                        vstar = bellman
                        astar = iia
                    end
                end
                vnew[ie, ia, id] = vstar
                gidx.a[ie, ia, id] = astar
                gidx.d[ie, ia, id] = id  # Moved inside to avoid extra loop
            end
        end
    end

    return vnew, gidx
end
