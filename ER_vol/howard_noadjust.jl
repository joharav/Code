function howard_noadjust(queuelong::Array{Float64}, util::Array{Float64}, iid, gidx::dtp.Ipol)
    vnew = zeros(sz.ne, sz.na, sz.nd)
    beta = pea[1]
    iiid = zeros(Int, sz.nd)

    Threads.@threads for id in 1:sz.nd
        iiid = iid[id]
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                iia = gidx.a[ie, ia, id]
                vnew[ie, ia, id] = util[ie, ia, id, iia] + beta * queuelong[ie, iia, iiid]
            end
        end
    end

    return vnew, gidx
end