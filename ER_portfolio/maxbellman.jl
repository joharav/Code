function maxbellman(queuelong::Array{Float64}, util::Array{Float64}, beta::Float64)
    vnew = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    gidx = dtp.Ipol(
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)), # a′  (foreign)
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)), # aa′ (local)
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd))  # d′
    )

    @Threads.threads for id in 1:sz.nd
        @Threads.threads for ia in 1:sz.na
            @Threads.threads for iaa in 1:sz.na
                @Threads.threads for iy in 1:sz.ny
                    @Threads.threads for ie in 1:sz.ne
                        vstar = -Inf; astar = 1; aastar = 1; dstar = 1
                        for iiaa in 1:sz.npa
                            for iia in 1:sz.npa
                                for iid in 1:sz.npd
                                    bellman = util[ie,iy,iaa,ia,id,iiaa,iia,iid] +
                                              beta*queuelong[ie,iy,iiaa,iia,iid]
                                    if bellman > vstar
                                        vstar = bellman; astar = iia; aastar = iiaa; dstar = iid
                                    end
                                end
                            end
                        end
                        gidx.a[ie,iy,iaa,ia,id] = astar
                        gidx.aa[ie,iy,iaa,ia,id] = aastar
                        gidx.d[ie,iy,iaa,ia,id] = dstar
                        vnew[ie,iy,iaa,ia,id]    = vstar
                    end
                end
            end
        end
    end
    return vnew, gidx
end
