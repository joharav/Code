function howard_noadjust(queuelong::Array{Float64,5},
    util::Array{Float64,8},
    iid_map::Vector{Int},
    old_iidx::dtp.Ipol, beta::Float64)
# queuelong: (ne,ny,npa,npa,npd)
# util:      (ne,ny,na,na,nd,npa,npa,npd) but only the slice at iid_map[id] is valid
    vnew = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)

    @Threads.threads for ie in 1:sz.ne
        for iy in 1:sz.ny, iaa in 1:sz.na, ia in 1:sz.na, id in 1:sz.nd
            iiaa = old_iidx.aa[ie,iy,iaa,ia,id]      # 1:sz.npa
            iia  = old_iidx.a[ie,iy,iaa,ia,id]       # 1:sz.npa
            iid  = iid_map[id]                       # pinned durable index

            vnew[ie,iy,iaa,ia,id] =
            util[ie,iy,iaa,ia,id,iiaa,iia,iid] + beta*queuelong[ie,iy,iiaa,iia,iid]
        end
    end

    # keep the durable policy index pinned
     old_iidx.d[:, :, :, :, :] .= reshape(iid_map, (1,1,1,1,sz.nd))
    return vnew, old_iidx
end
