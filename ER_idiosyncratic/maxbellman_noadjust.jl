function maxbellman_noadjust(queuelong::Array{Float64},util::Array{Float64},iid::Vector{Int64})
    beta        = pea[1]

    vnew = zeros(sz.ne,sz.ny,sz.na,sz.nd);
    gidx = dtp.Ipol(Int.(zeros(sz.ne,sz.ny,sz.na,sz.nd)), Int.(zeros(sz.ne,sz.ny,sz.na,sz.nd)))
    
    @Threads.threads for id in 1:sz.nd
        @Threads.threads for ia in 1:sz.na;
            @Threads.threads for iy in 1:sz.ny;
                @Threads.threads for ie in 1:sz.ne;
                    astar = 0; 
                    vstar = -Inf;
                    idd = iid[id]
                    for iia in 1:sz.npa;
                            bellman = util[ie,iy,ia,id,iia,idd] + beta*queuelong[ie,iy,iia,idd]
                            if bellman > vstar
                                vstar = bellman
                                astar = iia
                            end
                            if isnan(vstar)
                                println("NaN found in non-adjust Bellman!")
                            end
                            if bellman < -1e10
                                println("Suspiciously low Bellman: ", bellman, " at ie=", ie, " iy=", iy, " ia=", ia, " id=", id)
                            end
                            

                    end
                    gidx.a[ie,iy,ia,id] = astar
                    vnew[ie,iy,ia,id] = vstar
                end
            end
        end
    end


    if maximum(gidx.a) > sz.npa
        println("aint out of bounds")
    end

    return vnew::Array{Float64},gidx::dtp.Ipol
end