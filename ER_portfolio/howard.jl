function howard(queuelong::Array{Float64,5}, util::Array{Float64,8}, old_iidx::dtp.Ipol, beta::Float64)
    # queuelong: (ne,ny,npa,npa,npd)
    # util:      (ne,ny,na,na,nd,npa,npa,npd)  at the chosen (iiaa,iia,iid)
    vnew = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)

    @Threads.threads for ie in 1:sz.ne
        for iy in 1:sz.ny, iaa in 1:sz.na, ia in 1:sz.na, id in 1:sz.nd
            iiaa = old_iidx.aa[ie,iy,iaa,ia,id]   # ∈ 1:sz.npa
            iia  = old_iidx.a[ie,iy,iaa,ia,id]    # ∈ 1:sz.npa
            iid  = old_iidx.d[ie,iy,iaa,ia,id]    # ∈ 1:sz.npd

            @inbounds vnew[ie,iy,iaa,ia,id] =
                util[ie,iy,iaa,ia,id,iiaa,iia,iid] + beta*queuelong[ie,iy,iiaa,iia,iid]
        end
    end
    return vnew, old_iidx  # evaluation only; no policy change here
end
