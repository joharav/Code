function tinybellman(q::Array{Float64}, pr::Array{Float64}, old_gidx::dtp.Ipol)
    beta = pea[1]

    global sapol = zeros(Int, sz.nz, sz.ne, sz.na, sz.nd)
    global sdpol = zeros(Int, sz.nz, sz.ne, sz.na, sz.nd)
    global vnew = zeros(Float64, sz.nz, sz.ne, sz.na, sz.nd)
    gidx = deepcopy(old_gidx)

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for iz in 1:sz.nz
                    ja = old_gidx.a[iz, ie, ia, id]
                    jd = old_gidx.d[iz, ie, ia, id]
                    u_a = ja + sz.pad + max(1 + sz.pad - ja, 0) - max(ja + sz.pad - sz.npa, 0)
                    l_a = ja - sz.pad + max(1 + sz.pad - ja, 0) - max(ja + sz.pad - sz.npa, 0)
                    u_d = jd + sz.pad + max(1 + sz.pad - jd, 0) - max(jd + sz.pad - sz.npd, 0)
                    l_d = jd - sz.pad + max(1 + sz.pad - jd, 0) - max(jd + sz.pad - sz.npd, 0)
                    avsmall = beta * q[iz, ie, l_a:u_a, l_d:u_d] + pr[iz, ie, ia, id, l_a:u_a, l_d:u_d]
                    tmp = argmax(avsmall)
                    iia = tmp[1]
                    iid = tmp[2]
                    sapol[iz, ie, ia, id] = iia
                    sdpol[iz, ie, ia, id] = iid
                    vnew[iz, ie, ia, id] = avsmall[iia, iid]
                end
            end
        end
    end

    gidx.a = Int.(old_gidx.a) .+ Int.(sapol) .- sz.pad .- 1 .+ max.(sz.pad .+ 1 .- Int.(old_gidx.a), 0) .- max.(Int.(old_gidx.a) .+ sz.pad .- sz.npa, 0)
    gidx.d = Int.(old_gidx.d) .+ Int.(sdpol) .- sz.pad .- 1 .+ max.(sz.pad .+ 1 .- Int.(old_gidx.d), 0) .- max.(Int.(old_gidx.d) .+ sz.pad .- sz.npd, 0)

    return vnew::Array{Float64}, gidx::dtp.Ipol
end