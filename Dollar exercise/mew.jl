function makemew(answ::NamedTuple)
    k_wgt_a = inbetween(answ.g, true, :a) # True is for log spacing
    k_wgt_d = inbetween(answ.g, true, :d) # True is for log spacing

    interp_wgt, coarse_idx = makeweight(answ, k_wgt_a, k_wgt_d)

    mew = ones(sz.np, sz.ne, sz.nw, sz.na, sz.nd) ./ float(sz.np * sz.ne * sz.nw * sz.na * sz.nd)
    mew_new = mew

    for id in 1:sz.maxditer
        mew_new = makedist(coarse_idx, answ.i, interp_wgt, answ.g, mew_new, sz)
        disterr = maximum(abs.(mew .- mew_new))
        if disterr < sz.distol
            break
        else
            mew = mew_new
        end
    end
    mew = mew ./ sum(mew) # Rounding error
    return mew::Array{Float64, 5}
end

function makedist(coarse_idx, gidx, interp_wgt, grids, mew, sz)
    tmat = grids.t

    # Initialize the new distribution
    mew_old = mew

    for ip in 1:sz.np
        for ie in 1:sz.ne
            for iw in 1:sz.nw
                for ia in 1:sz.na
                    for id in 1:sz.nd
                        if mew_old[ip, ie, iw, ia, id] > 0.0
                            a_weight = interp_wgt[ip, ie, iw, ia, id].a_weight
                            d_weight = interp_wgt[ip, ie, iw, ia, id].d_weight
                            aidx = interp_wgt[ip, ie, iw, ia, id].a_lo
                            didx = interp_wgt[ip, ie, iw, ia, id].d_lo
                            a_increm = min(sz.na - aidx, 1)
                            d_increm = min(sz.nd - didx, 1)

                            for jp in 1:sz.np
                                for je in 1:sz.ne
                                    for jw in 1:sz.nw
                                        if tmat[jp, ip, je, ie, jw, iw] > 0.0
                                            mew[jp, je, jw, aidx, didx] += tmat[jp, ip, je, ie, jw, iw] * mew_old[ip, ie, iw, ia, id] * (1.0 - a_weight) * (1.0 - d_weight)
                                            if sz.na > sz.na || sz.nd > sz.nd
                                                mew[jp, je, jw, aidx + a_increm, didx] += tmat[jp, ip, je, ie, jw, iw] * mew_old[ip, ie, iw, ia, id] * a_weight * (1.0 - d_weight)
                                                mew[jp, je, jw, aidx, didx + d_increm] += tmat[jp, ip, je, ie, jw, iw] * mew_old[ip, ie, iw, ia, id] * (1.0 - a_weight) * d_weight
                                                mew[jp, je, jw, aidx + a_increm, didx + d_increm] += tmat[jp, ip, je, ie, jw, iw] * mew_old[ip, ie, iw, ia, id] * a_weight * d_weight
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return mew
end

function makeweight(answ::NamedTuple, k_wgt_a::NamedTuple, k_wgt_d::NamedTuple) 
    gidx  = answ.i
    pol   = answ.k 
    grids = answ.g

    kpg_a = grids.ap
    kpg_d = grids.dp

    interp_wgt = Dict()
    coarse_idx = Dict()

    for ip in 1:sz.np
        for ie in 1:sz.ne
            for iw in 1:sz.nw
                for ia in 1:sz.na
                    for id in 1:sz.nd
                        interp_wgt[ip, ie, iw, ia, id] = (a_weight = k_wgt_a.w[ia], d_weight = k_wgt_d.w[id], a_lo = k_wgt_a.lo[ia], d_lo = k_wgt_d.lo[id])
                        coarse_idx[ip, ie, iw, ia, id] = (a_idx = k_wgt_a.lo[ia], d_idx = k_wgt_d.lo[id])
                    end
                end
            end
        end
    end

    return interp_wgt, coarse_idx
end