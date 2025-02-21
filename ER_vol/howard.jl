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
    vnew = zeros(sz.ne,sz.na,sz.nd);
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na;
            Threads.@threads for ie in 1:sz.ne
                iia = old_gidx.a[ie,ia,id]; 
                iid = old_gidx.d[ie,ia,id];
                vnew[ie,ia,id] = beta*queuelong[ie,iia,iid] + util[ie,ia,id,iia,iid]   
            end
        end 
    end
    return vnew::Array{Float64}, old_gidx::dtp.Ipol; 
end
