function maxbellman(queuelong::Array{Float64},util::Array{Float64})
    vnew = zeros(sz.np,sz.ne,sz.nw,sz.na,sz.nd);

    gidx = dtp.Ipol( Int.(zeros(sz.np,sz.ne,sz.nw,sz.na,sz.nd)), Int.(zeros(sz.np,sz.ne,sz.nw,sz.na,sz.nd)))
    beta = pea[1]

    @Threads.threads for id in 1:sz.nd
        @Threads.threads for ia in 1:sz.na;
            @Threads.threads for iw in 1:sz.nw;
                @Threads.threads for ie in 1:sz.ne;
                    @Threads.threads for ip in 1:sz.np;
                        dstar = 0;
                        astar = 0; 
                        vstar = -Inf;
                        for iia in 1:sz.npa;
                            for iid in 1:sz.npd;
                                bellman = util[ip,ie,iw,ia,id,iia,iid] + beta*queuelong[ip,ie,iw,iia,iid]
                                if bellman > vstar
                                    vstar = bellman
                                    astar = iia
                                    dstar = iid
                                end
                            end
                        end
                        gidx.a[ip,ie,iw,ia,id] = astar
                        gidx.d[ip,ie,iw,ia,id] = dstar
                        vnew[ip,ie,iw,ia,id] = vstar
                    end
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