function tinybellman(q::Array{Float64,5}, pr::Array{Float64,8}, old_gidx::dtp.Ipol)
    # q  :: (ne, ny, iiaa, iia, iid)  expected value over z' at policy grids
    # pr :: (ne, ny, iaa,  ia,  id,   iiaa, iia, iid) current-period utility
    @assert size(q)  == (sz.ne, sz.ny, sz.npa, sz.npa, sz.npd)
    @assert size(pr) == (sz.ne, sz.ny, sz.na,  sz.na,  sz.nd,  sz.npa, sz.npa, sz.npd)

    beta = pea[1]

    vnew = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    gidx = deepcopy(old_gidx)

    # Flatten outer loops to avoid nested @threads
    total = sz.ne * sz.ny * sz.na * sz.na * sz.nd
    Threads.@threads for J in 1:total
        # Decode J -> (ie, iy, iaa, ia, id)
        t = J - 1
        id  = 1 + (t % sz.nd);  t ÷= sz.nd
        ia  = 1 + (t % sz.na);  t ÷= sz.na
        iaa = 1 + (t % sz.na);  t ÷= sz.na
        iy  = 1 + (t % sz.ny);  t ÷= sz.ny
        ie  = 1 +  t

        # previous policy indices
        ja  = gidx.a[ie,iy,iaa,ia,id]
        jaa = gidx.aa[ie,iy,iaa,ia,id]
        jd  = gidx.d[ie,iy,iaa,ia,id]

        # local windows (clamped)
        l_a  = max(ja  - sz.pad, 1);  u_a  = min(ja  + sz.pad, sz.npa)
        l_aa = max(jaa - sz.pad, 1);  u_aa = min(jaa + sz.pad, sz.npa)
        l_d  = max(jd  - sz.pad, 1);  u_d  = min(jd  + sz.pad, sz.npd)

        best = -Inf
        biia  = l_a
        biiaa = l_aa
        biid  = l_d

        @inbounds for iid in l_d:u_d
            # expected value term depends only on (iiaa,iia,iid)
            for iia in l_a:u_a
                for iiaa in l_aa:u_aa
                    val = pr[ie,iy,iaa,ia,id, iiaa,iia,iid] + beta * q[ie,iy, iiaa,iia,iid]
                    if val > best
                        best  = val
                        biia  = iia
                        biiaa = iiaa
                        biid  = iid
                    end
                end
            end
        end

        vnew[ie,iy,iaa,ia,id] = best
        gidx.a[ie,iy,iaa,ia,id] = biia
        gidx.aa[ie,iy,iaa,ia,id] = biiaa
        gidx.d[ie,iy,iaa,ia,id] = biid
    end

    return vnew::Array{Float64,5}, gidx::dtp.Ipol
end
