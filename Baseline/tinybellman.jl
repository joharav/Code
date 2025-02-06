function tinybellman(q::Array{Float64}, pr::Array{Float64}, old_gidx::dtp.Ipol)
    # This function maximizes the Bellman equation over a limited part of the policy space. 
    # Inputs

    # q: Array{Float64} - the value function
    # pr: Array{Float64} - the return function
    # old_gidx: Ipol - the old policy function

    # Outputs 

    # vnew: Array{Float64} - the new value function
    # gidx: Ipol - the new policy function

    beta = pea[1]
    global sapol = zeros(Int, sz.np, sz.na, sz.nd)
    global sdpol = zeros(Int, sz.np, sz.na, sz.nd)
    global  vnew = zeros(Float64, sz.np, sz.na, sz.nd)

    gidx = deepcopy(old_gidx)

    # NOTE!!! This will not work if the last best point is too close to the end of the state space. 
    # It will throw an error. 
    # This is good. It means you have a bad state space. 
    for id in 1:sz.nd
        for ia in 1:sz.na
           for ip in 1:sz.np
               ja = old_gidx.a[ip,ia,id] 
               jd = old_gidx.d[ip,ia,id] 
               u_a = ja + sz.pad + max(1+sz.pad-ja,0) - max(ja+sz.pad-sz.npa, 0);
               l_a = ja - sz.pad + max(1+sz.pad-ja,0) - max(ja+sz.pad-sz.npa, 0);
               u_d = jd + sz.pad + max(1+sz.pad-jd,0) - max(jd+sz.pad-sz.npd, 0);
               l_d = jd - sz.pad + max(1+sz.pad-jd,0) - max(jd+sz.pad-sz.npd, 0);
               avsmall = beta * q[ip, l_a:u_a, l_d:u_d] + pr[ip,ia,id,l_a:u_a, l_d:u_d] 
               tmp = argmax(avsmall);
               iia = tmp[1];
               iid = tmp[2]; 
               sapol[ip,ia,id] = iia
               sdpol[ip,ia,id] = iid 
               vnew[ip,ia,id]  = avsmall[iia,iid]
           end 
        end  
    end 

    gidx.a = Int.(old_gidx.a) .+ Int.(sapol) .- sz.pad .- 1 .+ max.(sz.pad .+ 1 .- Int.(old_gidx.a),0) .- max.(Int.(old_gidx.a) .+ sz.pad .- sz.npa,0)
    gidx.d = Int.(old_gidx.d) .+ Int.(sdpol) .- sz.pad .- 1 .+ max.(sz.pad .+ 1 .- Int.(old_gidx.d),0) .- max.(Int.(old_gidx.d) .+ sz.pad .- sz.npd,0)

    return vnew::Array{Float64}, gidx::dtp.Ipol
end 