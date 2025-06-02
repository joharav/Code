function maxbellman_noadjust(queuelong::Array{Float64},util::Array{Float64},iid::Vector{Int64})
    beta        = pea[1]

    vnew = zeros(sz.ne,sz.na,sz.nd);
    gidx = dtp.Ipol(Int.(zeros(sz.ne,sz.na,sz.nd)), Int.(zeros(sz.ne,sz.na,sz.nd)))
    
    @Threads.threads for id in 1:sz.nd
        @Threads.threads for ia in 1:sz.na;
            @Threads.threads for ie in 1:sz.ne;
                astar = 0; 
                vstar = -Inf;
                idd = iid[id]
                for iia in 1:sz.npa;
                        bellman = util[ie,ia,id,iia,idd] + beta*queuelong[ie,iia,idd]
                        if bellman > vstar
                            vstar = bellman
                            astar = iia
                        end
                end
                gidx.a[ie,ia,id] = astar
                vnew[ie,ia,id] = vstar
            end
        end
    end


    if maximum(gidx.a) > sz.npa
        println("aint out of bounds")
    end

    return vnew::Array{Float64},gidx::dtp.Ipol
end