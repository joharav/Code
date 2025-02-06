function makemew(answ::NamedTuple)

    a_wgt = inbetween(answ.g, true, :a)  # Assuming answ.g.k is the grid for a and true indicates log spacing
    d_wgt = inbetween(answ.g, true, :d)  # Assuming answ.g.d is the grid for d and true indicates log spacing

    interp_wgt, coarse_idx = makeweight(answ, a_wgt, d_wgt)

    mew = ones(sz.np, sz.na, sz.nd) ./ float(sz.np * sz.na * sz.nd)
    mew_new = mew

    for id in 1:sz.maxditer
        mew_new = makedist(coarse_idx, answ.gidx, interp_wgt, answ.g, mew_new)
        disterr = maximum(abs.(mew .- mew_new))
        if disterr < sz.distol
            break
        else
            mew = mew_new
        end
    end

    mew = mew ./ sum(mew)  # Rounding error
    return mew::Array{Float64, 3}  # Return a 3D array (matrix) now
end

function makedist(coarse_idx, gidx::dtp.Ipol, interp_wgt::Array{Float64, 5}, grids::Vector{Float64}, mew::Array{Float64, 3})

    tmat = grids.t

    # Initialize the new distribution
    mew_old = mew

    new_mew = zeros(size(mew))

    for ia in 1:sz.na
        for ip in 1:sz.np
            for id in 1:sz.nd
                if mew_old[ip, ia, id] > 0.0
                    a_weight = interp_wgt[ip, ia, id, Int(gidx.a[ip, ia, id]), Int(gidx.d[ip, ia, id])]
                    a_increm = min(sz.na - Int(coarse_idx.a[ip, ia, id]), 1)
                    d_weight = interp_wgt[ip, ia, id, Int(gidx.a[ip, ia, id]), Int(gidx.d[ip, ia, id])]
                    d_increm = min(sz.nd - Int(coarse_idx.d[ip, ia, id]), 1)
                    for ipp in 1:sz.np
                        if tmat[ip, ipp] > 0.0
                            aidx = Int(coarse_idx.a[ip, ia, id])
                            didx = Int(coarse_idx.d[ip, ia, id])
                            new_mew[ipp, aidx, didx] += tmat[ipp, ip] * mew_old[ip, ia, id] * a_weight * d_weight
                            if sz.npa > sz.na && sz.npd > sz.nd
                                new_mew[ipp, aidx + a_increm, didx + d_increm] += tmat[ipp, ip] * mew_old[ip, ia, id] * (1.0 - a_weight) * (1.0 - d_weight)
                            end
                        end
                    end
                end
            end
        end
    end

    return new_mew::Array{Float64, 3}
end


function makeweight(answ::NamedTuple, a_wgt::NamedTuple, d_wgt::NamedTuple) 
    gidx  = answ.gidx
    pol   = answ.pol
    grids = answ.g  # Changed from answ.grids to answ.g

    apg   = grids.ap
    dpg   = grids.dp

    interp_wgt = zeros(sz.np, sz.na, sz.nd, sz.npa, sz.npd)
    coarse_idx = (a=Int.(zeros(sz.np, sz.na, sz.nd)), d=Int.(zeros(sz.np, sz.na, sz.nd)))

    if (sz.na == sz.npa && sz.nd == sz.npd)  # Changed sz.npk to sz.npa for consistency
        for ia in 1:sz.na
            for ip in 1:sz.np
                for id in 1:sz.nd
                    coarse_idx.a[ip, ia, id] = gidx.a[ip, ia, id]  # Changed gidx[ip, ia, id] to gidx.a[ip, ia, id]
                    coarse_idx.d[ip, ia, id] = gidx.d[ip, ia, id]  # Changed gidx[ip, ia, id] to gidx.d[ip, ia, id]

                    for iia in 1:sz.npa  # Changed sz.npk to sz.npa
                        for iid in 1:sz.npd
                            if (pol.a[ip, ia, id] * pol.d[ip, ia, id]) == (apg[iia] * dpg[iid])
                                interp_wgt[ip, ia, id, iia, iid] = 1.0
                            end
                        end
                    end
                end
            end
        end
    else
        for ia in 1:sz.na
            for ip in 1:sz.np
                for id in 1:sz.nd
                    coarse_idx.a[ip, ia, id] = a_wgt.lo[gidx.a[ip, ia, id]]  # Changed gidx[ip, ia, id] to gidx.a[ip, ia, id]
                    coarse_idx.d[ip, ia, id] = d_wgt.lo[gidx.d[ip, ia, id]]  # Changed gidx[ip, ia, id] to gidx.d[ip, ia, id]
                end
            end
        end
        for iia in 1:sz.npa
            for iid in 1:sz.npd
                for ia in 1:sz.na
                    for ip in 1:sz.np
                        for id in 1:sz.nd
                            if (pol.a[ip, ia, id] * pol.d[ip, ia, id]) == (apg[iia] * dpg[iid])
                                interp_wgt[ip, ia, id, iia, iid] = a_wgt.w[iia] * d_wgt.w[iid]
                            end
                        end
                    end
                end
            end
        end
    end

    return interp_wgt::Array{Float64, 5}, coarse_idx
end
