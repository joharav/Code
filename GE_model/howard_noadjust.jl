function howard_noadjust(queuelong::Array{Float64}, util::Array{Float64}, iid, gidx::dtp.Ipol)
    vnew = zeros(sz.nz, sz.ne, sz.na, sz.nd)
    beta = pea[1]
    iiid = zeros(Int, sz.nd)

    Threads.@threads for id in 1:sz.nd
        iiid = iid[id]
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for iz in 1:sz.nz
                    iia = gidx.a[iz, ie, ia, id]
                    vnew[iz, ie, ia, id] = util[iz, ie, ia, id, iia,iiid] + beta * queuelong[iz, ie, iia, iiid]
                end
            end
        end
    end

    return vnew, gidx
end