function howard_noadjust(queuelong::Array{Float64}, util::Array{Float64}, gidx::dtp.Ipol)
    vnew = zeros(sz.ne, sz.na, sz.nd)
    beta = pea[1]

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                iia = gidx.a[ie, ia, id]
                vnew[ie, ia, id] = util[ie, ia, id, iia] + beta * queuelong[ie, iia, id]
                gidx.d[ie, ia, id] = id
            end
        end
    end

    return vnew, gidx
end