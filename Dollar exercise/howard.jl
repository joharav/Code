function howard(queuelong::Array{Float64},util::Array{Float64},old_iidx::dtp.Ipol);

    # Inputs
    # queuelong: interpolated value function
    # util: utility array
    # old_aint current policy integer function
    # old_dint current policy integer function

    # Outputs
    # vnew: new value function
    beta = pea[1]
    old_gidx = deepcopy(old_iidx)
    vnew = zeros(sz.np,sz.ne,sz.nw,sz.na,sz.nd);
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na;
            Threads.@threads for iw in 1:sz.nw
                Threads.@threads for ie in 1:sz.ne
                    Threads.@threads for ip in 1:sz.np;
                        iia = old_iidx.a[ip,ie,iw,ia,id]; 
                        iid = old_iidx.d[ip,ie,iw,ia,id];
                        vnew[ip,ie,iw,ia,id] = beta*queuelong[ip,ie,iw,iia,iid] + util[ip,ie,iw,ia,id,iia,iid]   
                    end
                end
            end
        end 
    end
    return vnew::Array{Float64}, old_gidx::dtp.Ipol; 
end
