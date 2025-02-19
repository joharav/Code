function makemew(answ::NamedTuple)
    k_wgt_a = inbetween(answ.g, true, :a) # True is for log spacing
    k_wgt_d = inbetween(answ.g, true, :d) # True is for log spacing

    interp_wgt, coarse_idx = makeweight(answ, k_wgt_a, k_wgt_d)

    mew = ones(sz.np, sz.ne, sz.na, sz.nd) ./ float(sz.np * sz.ne * sz.na * sz.nd)
    mew_new = mew

    for id in 1:sz.maxditer
        mew_new = makedist(coarse_idx, interp_wgt, answ.g, mew_new)
        disterr = maximum(abs.(mew .- mew_new))
        if disterr < sz.distol
            break
        else
            mew = mew_new
        end
    end
    mew = mew ./ sum(mew) # Rounding error
    return mew::Array{Float64, 4}
end

function makedist(coarse_idx::Dict{Any, Any}, interp_wgt::Dict{Any, Any}, grids::NamedTuple, mew::Array{Float64, 4})
    tmat = grids.t
    tmat_e = tmat[1:sz.ne, 1:sz.ne]
    tmat_p = tmat[1:sz.ne:end, 1:sz.ne:end]

    # Initialize the new distribution
    mew_old = mew
    mew_new = zeros(size(mew))

    # Loop over all state variables to update the distribution
    for id in 1:sz.nd
        for ia in 1:sz.na
            for ie in 1:sz.ne
                for ip in 1:sz.np
                    # Get the coarse indices and interpolation weights
                    (lo_a, hi_a, lo_d, hi_d) = coarse_idx[(ip, ie, ia, id)]
                    (w_a, w_d) = interp_wgt[(ip, ie, ia, id)]

                    # Update the distribution
                    for a_increm in 0:1
                        for d_increm in 0:1
                            a_weight = w_a[a_increm + 1]
                            d_weight = w_d[d_increm + 1]
                            aidx = lo_a + a_increm
                            didx = lo_d + d_increm

                            # Ensure indices are within bounds
                            if aidx <= sz.na && didx <= sz.nd
                                for jp in 1:sz.np
                                    for je in 1:sz.ne
                                        tmat_val = tmat_e[ie, je] * tmat_p[ip, jp]
                                        mew_new[jp, je, aidx, didx] += tmat_val * mew_old[ip, ie, ia, id] * a_weight * d_weight
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return mew_new
end

function makeweight(answ::NamedTuple, k_wgt_a::NamedTuple, k_wgt_d::NamedTuple)
    # Extract data components from the named tuples
    gidx  = answ.gidx
    pol   = answ.pol
    grids = answ.g
    kpg_a = grids.ap  # Grid points for `a`
    kpg_d = grids.dp  # Grid points for `d` policy 

    # Initialize interpolation weights and coarse indices dictionaries
    interp_wgt = Dict()
    coarse_idx = Dict()

    # Loop over all state variables to calculate weights and indices
    for id in 1:sz.nd
        for ia in 1:sz.na
            for ie in 1:sz.ne
                for ip in 1:sz.np
                    # Get the policy indices
                    a_idx = gidx.a[ip, ie, ia, id]
                    d_idx = gidx.d[ip, ie, ia, id]

                    # Calculate interpolation weights for assets
                    lo_a = k_wgt_a.lo[a_idx]
                    hi_a = k_wgt_a.hi[a_idx]
                    w_a = [k_wgt_a.w[a_idx], 1.0 - k_wgt_a.w[a_idx]]

                    # Calculate interpolation weights for durables
                    lo_d = k_wgt_d.lo[d_idx]
                    hi_d = k_wgt_d.hi[d_idx]
                    w_d = [k_wgt_d.w[d_idx], 1.0 - k_wgt_d.w[d_idx]]

                    # Ensure indices are within bounds
                    lo_a = min(max(lo_a, 1), sz.na)
                    hi_a = min(max(hi_a, 1), sz.na)
                    lo_d = min(max(lo_d, 1), sz.nd)
                    hi_d = min(max(hi_d, 1), sz.nd)

                    # Store the weights and indices
                    interp_wgt[(ip, ie, ia, id)] = (w_a, w_d)
                    coarse_idx[(ip, ie, ia, id)] = (lo_a, hi_a, lo_d, hi_d)
                end
            end
        end
    end

    return interp_wgt, coarse_idx
end