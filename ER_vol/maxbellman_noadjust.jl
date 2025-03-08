function maxbellman_noadjust(queuelong::Array{Float64}, util::Array{Float64},iid)
    vnew = zeros(sz.ne,sz.na,sz.nd);
    gidx = dtp.Ipol(Int.(zeros(sz.ne,sz.na,sz.nd)), Int.(zeros(sz.ne,sz.na,sz.nd)))
    beta = pea[1]

    Threads.@threads for id in 1:sz.nd
        for ia in 1:sz.na
            for ie in 1:sz.ne
                dstar = 0;
                astar = 0; 
                vstar = -Inf;
                # for iia in 1:sz.npa
                #     bellman = util[ie, ia, id, iia] + beta * queuelong[ie, iia, iiid]
                #     if bellman > vstar
                #         vstar = bellman
                #         astar = iia
                #         dstar = iiid
                #     end
                # end
                for iia in 1:sz.npa;
                    for iid in 1:sz.npd;
                        bellman = util[ie,ia,id,iia,iid] + beta*queuelong[ie,iia,iid]
                        if bellman > vstar
                            vstar = bellman
                            astar = iia
                            dstar = iid
                        end
                    end
                end


                vnew[ie, ia, id] = vstar
                gidx.a[ie, ia, id] = astar
                gidx.d[ie, ia, id] = dstar  # Moved inside to avoid extra loop
            end
        end
    end
    if maximum(gidx.a) > sz.npa
        println("aint out of bounds")
    end
    if maximum(gidx.d) > sz.npd
        println("dint out of bounds")
    end
    return vnew, gidx
end
