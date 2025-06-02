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
    vnew = zeros(sz.nz,sz.ne,sz.na,sz.nd);
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na;
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for iz in 1:sz.nz
                    iia = old_gidx.a[iz,ie,ia,id]; 
                    iid = old_gidx.d[iz,ie,ia,id];
                    vnew[iz,ie,ia,id] = beta*queuelong[iz,ie,iia,iid] + util[iz,ie,ia,id,iia,iid]   
                end
            end
        end 
    end
    return vnew::Array{Float64}, old_gidx::dtp.Ipol; 
end
