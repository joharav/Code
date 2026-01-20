# ==========================================================================
# 4D MODEL: Fill in value function from state grids to policy grids
# obj:     (ne, ny, nw, nd)
# objlong: (ne, ny, npw, npd)
# ==========================================================================

function fillin(obj::Array{Float64,4}, g::NamedTuple)
    objlong = Array{Float64}(undef, sz.ne, sz.ny, sz.npw, sz.npd)

    # Precompute brackets for policy grids on state grids (fast + correct)
    wL = Vector{Int}(undef, sz.npw)
    wU = Vector{Int}(undef, sz.npw)
    ww = Vector{Float64}(undef, sz.npw)
    @inbounds for iwp in 1:sz.npw
        l,u,w = brack1d(g.w, g.wp[iwp])
        wL[iwp] = l; wU[iwp] = u; ww[iwp] = w
    end

    dL = Vector{Int}(undef, sz.npd)
    dU = Vector{Int}(undef, sz.npd)
    wd = Vector{Float64}(undef, sz.npd)
    @inbounds for idp in 1:sz.npd
        l,u,w = brack1d(g.d, g.dp[idp])
        dL[idp] = l; dU[idp] = u; wd[idp] = w
    end

    # Parallelize over (iwp,idp). Each thread writes disjoint slices => safe.
    Threads.@threads for J in 1:(sz.npw * sz.npd)
        iwp = 1 + (J - 1) % sz.npw
        idp = 1 + (J - 1) ÷ sz.npw

        iwL = wL[iwp]; iwU = wU[iwp]; wwt = ww[iwp]
        idL = dL[idp]; idU = dU[idp]; dwt = wd[idp]

        @inbounds for iy in 1:sz.ny, ie in 1:sz.ne
            v00 = obj[ie, iy, iwL, idL]
            v10 = obj[ie, iy, iwU, idL]
            v01 = obj[ie, iy, iwL, idU]
            v11 = obj[ie, iy, iwU, idU]

            v0 = (1.0 - wwt) * v00 + wwt * v10
            v1 = (1.0 - wwt) * v01 + wwt * v11

            objlong[ie, iy, iwp, idp] = (1.0 - dwt) * v0 + dwt * v1
        end
    end

    return objlong
end
