function maxbellman_noadjust(queuelong::Array{Float64}, util::Array{Float64}, iid::Vector{Int64})
    beta = pea[1]
    vnew = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    gidx = dtp.Ipol(
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)),
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)),
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd))
    )

    @Threads.threads for id in 1:sz.nd
        @Threads.threads for ia in 1:sz.na
            @Threads.threads for iaa in 1:sz.na
                @Threads.threads for iy in 1:sz.ny
                    @Threads.threads for ie in 1:sz.ne
                        vstar = -Inf; astar = 1; aastar = 1
                        idd = iid[id]
                        for iiaa in 1:sz.npa
                            for iia in 1:sz.npa
                                bellman = util[ie,iy,iaa,ia,id,iiaa,iia,idd] +
                                          beta*queuelong[ie,iy,iiaa,iia,idd]
                                if bellman > vstar
                                    vstar = bellman; astar = iia; aastar = iiaa
                                end
                            end
                        end
                        gidx.a[ie,iy,iaa,ia,id] = astar
                        gidx.aa[ie,iy,iaa,ia,id] = aastar
                        gidx.d[ie,iy,iaa,ia,id] = idd
                        vnew[ie,iy,iaa,ia,id]    = vstar
                    end
                end
            end
        end
    end
    return vnew, gidx
end
