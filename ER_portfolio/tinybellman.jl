function tinybellman(q::Array{Float64}, pr::Array{Float64}, old_gidx::dtp.Ipol)
    beta = pea[1]

    sapol_a  = zeros(Int, sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    sapol_aa = zeros(Int, sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    sdpol    = zeros(Int, sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    vnew     = zeros(Float64, sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    gidx     = deepcopy(old_gidx)

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for iaa in 1:sz.na
                Threads.@threads for iy in 1:sz.ny
                    Threads.@threads for ie in 1:sz.ne
                        ja  = old_gidx.a[ie,iy,iaa,ia,id]
                        jaa = old_gidx.aa[ie,iy,iaa,ia,id]
                        jd  = old_gidx.d[ie,iy,iaa,ia,id]

                        u_a  = ja  + sz.pad + max(1 + sz.pad - ja,  0) - max(ja  + sz.pad - sz.npa, 0)
                        l_a  = ja  - sz.pad + max(1 + sz.pad - ja,  0) - max(ja  + sz.pad - sz.npa, 0)
                        u_aa = jaa + sz.pad + max(1 + sz.pad - jaa, 0) - max(jaa + sz.pad - sz.npa, 0)
                        l_aa = jaa - sz.pad + max(1 + sz.pad - jaa, 0) - max(jaa + sz.pad - sz.npa, 0)
                        u_d  = jd  + sz.pad + max(1 + sz.pad - jd,  0) - max(jd  + sz.pad - sz.npd, 0)
                        l_d  = jd  - sz.pad + max(1 + sz.pad - jd,  0) - max(jd  + sz.pad - sz.npd, 0)

                        avsmall = beta .* q[ie,iy,l_aa:u_aa, l_a:u_a, l_d:u_d] .+
                                  pr[ie,iy,iaa,ia, l_aa:u_aa, l_a:u_a, l_d:u_d]
                        tmp = argmax(avsmall)
                        iiaa = tmp[1]; iia = tmp[2]; iid = tmp[3]

                        sapol_aa[ie,iy,iaa,ia,id] = iiaa
                        sapol_a[ie,iy,iaa,ia,id] = iia
                        sdpol[ie,iy,iaa,ia,id] = iid
                        vnew[ie,iy,iaa,ia,id] = avsmall[iiaa,iia,iid]
                    end
                end
            end
        end
    end

    gidx.aa = Int.(old_gidx.aa) .+ Int.(sapol_aa) .- sz.pad .- 1 .+
              max.(sz.pad .+ 1 .- Int.(old_gidx.aa), 0) .- max.(Int.(old_gidx.aa) .+ sz.pad .- sz.npa, 0)
    gidx.a  = Int.(old_gidx.a ) .+ Int.(sapol_a ) .- sz.pad .- 1 .+
              max.(sz.pad .+ 1 .- Int.(old_gidx.a ), 0) .- max.(Int.(old_gidx.a ) .+ sz.pad .- sz.npa, 0)
    gidx.d  = Int.(old_gidx.d ) .+ Int.(sdpol   ) .- sz.pad .- 1 .+
              max.(sz.pad .+ 1 .- Int.(old_gidx.d ), 0) .- max.(Int.(old_gidx.d ) .+ sz.pad .- sz.npd, 0)

    return vnew, gidx
end
