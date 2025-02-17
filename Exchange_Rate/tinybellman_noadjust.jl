function tinybellman_noadjust(queuelong::Array{Float64}, util::Array{Float64}, gidx::dtp.Ipol)
    vnew = zeros(sz.np, sz.ne, sz.na, sz.nd)
    beta = pea[1]

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    iia = gidx.a[ip, ie, ia, id]
                    vnew[ip, ie, ia, id] = util[ip, ie, ia, id, iia] + beta * queuelong[ip, ie, iia, id]
                end
            end
        end
    end

    # Update gidx.d outside the loops
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    gidx.d[ip, ie, ia, id] = id
                end
            end
        end
    end

    return vnew, gidx
end