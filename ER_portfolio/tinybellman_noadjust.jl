function tinybellman_noadjust(q::Array{Float64,5}, pr::Array{Float64,8}, iid::Vector{Int}, old_gidx::dtp.Ipol)
    beta = pea[1]

    vnew = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    gidx = deepcopy(old_gidx)

    for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na, iy in 1:sz.ny, ie in 1:sz.ne
        ja  = old_gidx.a[ie,iy,iaa,ia,id]
        jaa = old_gidx.aa[ie,iy,iaa,ia,id]
        idd = iid[id]

        # local windows (clamped)
        u_a  = min(ja  + sz.pad, sz.npa); l_a  = max(ja  - sz.pad, 1)
        u_aa = min(jaa + sz.pad, sz.npa); l_aa = max(jaa - sz.pad, 1)

        # 2-D block over (iiaa,iia)
        best = -Inf
        best_iiaa = l_aa
        best_iia  = l_a
        for iiaa in l_aa:u_aa
            for iia in l_a:u_a
                val = beta * q[ie,iy,iiaa,iia,idd] + pr[ie,iy,iaa,ia,id,iiaa,iia,idd]
                if val > best
                    best      = val
                    best_iiaa = iiaa
                    best_iia  = iia
                end
            end
        end

        gidx.aa[ie,iy,iaa,ia,id] = best_iiaa
        gidx.a[ie,iy,iaa,ia,id] = best_iia
        # d' fixed already by iid
        gidx.d[ie,iy,iaa,ia,id] = idd
        vnew[ie,iy,iaa,ia,id]    = best
    end

    return vnew::Array{Float64,5}, gidx::dtp.Ipol
end
