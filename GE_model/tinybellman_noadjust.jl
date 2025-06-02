function tinybellman_noadjust(q::Array{Float64}, pr::Array{Float64}, iid::Vector{Int64}, old_gidx::dtp.Ipol)
    beta = pea[1]

    global sapol = zeros(Int, sz.nz, sz.ne, sz.na, sz.nd)
    global vnew = zeros(Float64, sz.nz, sz.ne, sz.na, sz.nd)
    gidx = deepcopy(old_gidx)

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for iz in 1:sz.nz
                    ja = old_gidx.a[iz, ie, ia, id]
                    idd = iid[id]

                    # Define bounds for asset policy search
                    u_a = ja + sz.pad + max(1 + sz.pad - ja, 0) - max(ja + sz.pad - sz.npa, 0)
                    l_a = ja - sz.pad + max(1 + sz.pad - ja, 0) - max(ja + sz.pad - sz.npa, 0)


                    # Value function maximization over a' only
                    avsmall = beta * q[iz, ie, l_a:u_a, idd] + pr[iz, ie, ia, id, l_a:u_a, idd]
                    tmp = argmax(avsmall)
                    iia = tmp[1]

                    # Store optimal choices
                    sapol[iz, ie, ia, id] = iia
                    vnew[iz, ie, ia, id] = avsmall[iia]
                end
            end
        end
    end

    # Update policy function: d' remains the same
    gidx.a = Int.(old_gidx.a) .+ Int.(sapol) .- sz.pad .- 1 .+ 
             max.(sz.pad .+ 1 .- Int.(old_gidx.a), 0) .- max.(Int.(old_gidx.a) .+ sz.pad .- sz.npa, 0)


    return vnew::Array{Float64}, gidx::dtp.Ipol
end