function maxbellman(queuelong::Array{Float64},util::Array{Float64})
    vnew = zeros(sz.ne,sz.ny,sz.na,sz.nd);
    gidx = dtp.Ipol(Int.(zeros(sz.ne,sz.ny,sz.na,sz.nd)), Int.(zeros(sz.ne,sz.ny,sz.na,sz.nd)))
    beta = pea[1]

    
    @Threads.threads for id in 1:sz.nd
        @Threads.threads for ia in 1:sz.na;
            @Threads.threads for iy in 1:sz.ny;
                @Threads.threads for ie in 1:sz.ne;
                    dstar = 0;
                    astar = 0; 
                    vstar = -Inf;
                    for iia in 1:sz.npa;
                        for iid in 1:sz.npd;
                            bellman = util[ie,iy,ia,id,iia,iid] + beta*queuelong[ie,iy,iia,iid]
                            if bellman > vstar
                                vstar = bellman
                                astar = iia
                                dstar = iid
                            end
                        end
                    end
                    gidx.a[ie,iy,ia,id] = astar
                    gidx.d[ie,iy,ia,id] = dstar
                    vnew[ie,iy,ia,id] = vstar
                end
            end
        end
    end
    if maximum(gidx.a) > sz.npa
        println("aint out of bounds")
    end
    if maximum(gidx.d) > sz.npd
        println("dint out of bounds")
    end
    return vnew::Array{Float64},gidx::dtp.Ipol
end