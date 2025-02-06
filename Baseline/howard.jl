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
    vnew = zeros(sz.np,sz.na,sz.nd);
    Threads.@threads for id in 1:sz.nd
        for ia in 1:sz.na;
             for ip in 1:sz.np;
                 iia = old_iidx.a[ip,ia,id]; 
                 iid = old_iidx.d[ip,ia,id];
                 vnew[ip,ia,id] = beta*queuelong[ip,iia,iid] + util[ip,ia,id,iia,iid]   
            end
        end 
    end
    return vnew::Array{Float64}, old_gidx::dtp.Ipol; 
end
