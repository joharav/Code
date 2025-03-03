# Define the function makemew, which takes a tuple as Inputs and returns a matrix as Output
function makemew(answ::NamedTuple)

    a_wgt = inbetween(answ.g,true,:a) # True is for log spacing
    d_wgt = inbetween(answ.g,true,:d) # True is for log spacing

    interp_wgt, coarse_idx = makeweight(answ, a_wgt, d_wgt)


    mew = ones(sz.ne,sz.na,sz.nd) ./ float(sz.ne*sz.na*sz.nd)
    mew_new = mew

    for id in 1:sz.maxditer
        mew_new = makedist(coarse_idx, interp_wgt, answ.g, mew_new)
        disterr = maximum(abs.(mew.-mew_new))
        if disterr < sz.distol
          break
        else
          mew = mew_new
        end
    end
    mew = mew./sum(mew) #Rounding error
    return mew::Matrix{Float64}

end

function makedist(coarse_idx, interp_wgt, grids, mew)

    tmat = grids.t
    mew_old = mew


    for id in 1:sz.nd
        for ia in 1:sz.na 
            for ie in 1:sz.ne
                if mew_old[ie,ia,id] > 0.0
                    a_weight = interp_wgt[ie,ia,id].a_weight
                    d_weight = interp_wgt[ie,ia,id].d_weight

                    aidx = Int(coarse_idx[ie,ia,id]).a_lo; 
                    didx = Int(coarse_idx[ie,ia,id]).d_lo; 

                    a_increm = min(sz.na-aidx,1)
                    d_increm = min(sz.nd-didx,1)



                    for iee in 1:sz.ne
                        if tmat[iee,ie] > 0.0
                            mew[iee,aidx,didx] = mew[iee,aidx,didx] + tmat[iee,ie]*mew_old[ie,ia,id] * (1-a_weight) * (1-d_weight)

                            if sz.npa > sz.na + sz.npd > sz.nd  
                                mew[iee,aidx+a_increm,didx] = mew[iee,aidx+a_increm,didx] + tmat[iee,ie]*mew_old[ie,ia,id] * a_weight * (1.0-d_weight)

                                mew[iee,aidx,didx+d_increm] = mew[iee,aidx,didx+d_increm] + tmat[iee,ie]*mew_old[ie,ia,id] * (1.0-a_weight) * d_weight

                                mew[iee,aidx+a_increm,didx+d_increm] = mew[iee,aidx+a_increm,didx+d_increm] + tmat[iee,ie]*mew_old[ie,ia,id] * a_weight * d_weight     
                            end                       
                        end
                    end
                end
            end
        end
    end
    
    return mew
end

function makeweight(answ::NamedTuple, a_wgt::NamedTuple, d_wgt::NamedTuple) 
    gidx  = answ.gidx


    # Initialize interpolation weight storage
    interp_wgt = [(a_weight=0.0, d_weight=0.0) for _ in 1:sz.ne, _ in 1:sz.na, _ in 1:sz.nd]
    
    # Initialize coarse index storage
    coarse_idx = [(a_lo=0, d_lo=0) for _ in 1:sz.ne, _ in 1:sz.na, _ in 1:sz.nd]

    for id in 1:sz.nd
        for ia in 1:sz.na
            for ie in 1:sz.ne
                g_idx = gidx[ie, ia, id]

                # Assign lower indices (grid positions just below policy choices)
                a_lo = a_wgt.lo[g_idx]
                d_lo = d_wgt.lo[g_idx]

                # Assign interpolation weights
                a_weight = a_wgt.w[g_idx]
                d_weight = d_wgt.w[g_idx]

                # Store in coarse_idx and interp_wgt
                coarse_idx[ie, ia, id] = (a_lo=a_lo, d_lo=d_lo)
                interp_wgt[ie, ia, id] = (a_weight=a_weight, d_weight=d_weight)
            end
        end
    end

    return interp_wgt, coarse_idx
end

