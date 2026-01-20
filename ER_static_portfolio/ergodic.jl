# ==========================================================================
# 4D MODEL: Compute ergodic distribution
# State: (e, y, w, d)
# ==========================================================================

function compute_ergodic(ans::NamedTuple; tol=sz.distol, max_iter=sz.maxditer)
    g = ans.g
    tmat = g.t  # size (ne*ny) × (ne*ny)
    
    ne, ny, nw, nd = sz.ne, sz.ny, sz.nw, sz.nd
    nzy = ne * ny
    nstates = ne * ny * nw * nd
    
    # Policies: 4D (ie, iy, iw, id)
    pol_w = ans.pol.w
    pol_d = ans.pol.d
    
    @assert size(pol_w) == (ne, ny, nw, nd)
    @assert size(pol_d) == (ne, ny, nw, nd)
    @assert size(tmat) == (nzy, nzy)
    
    # Linear index (column-major order)
    flat_idx(ie, iy, iw, id) = 
        ((((ie-1)*ny + (iy-1))*nw + (iw-1))*nd + (id-1)) + 1
    
    find_index(val, grid) = clamp(searchsortedfirst(grid, val), 1, length(grid))
    
    dist = fill(1.0 / nstates, nstates)
    dist_new = similar(dist)
    
    for it in 1:max_iter
        fill!(dist_new, 0.0)
        
        @inbounds for ie in 1:ne, iy in 1:ny, iw in 1:nw, id in 1:nd
            idx = flat_idx(ie, iy, iw, id)
            mass = dist[idx]
            mass == 0.0 && continue
            
            # Next-period policy indices
            iw′ = find_index(pol_w[ie, iy, iw, id], g.w)
            id′ = find_index(pol_d[ie, iy, iw, id], g.d)
            
            # Transition over z = (e, y)
            iz = (iy - 1) * ne + ie
            for iz′ in 1:nzy
                ie′ = mod1(iz′, ne)
                iy′ = div(iz′ - 1, ne) + 1
                prob = tmat[iz′, iz]
                prob == 0.0 && continue
                idx′ = flat_idx(ie′, iy′, iw′, id′)
                dist_new[idx′] += prob * mass
            end
        end
        
        err = maximum(abs.(dist_new .- dist))
        dist, dist_new = dist_new, dist
        if err < tol
            if settings.verbose
                println("Ergodic converged after $it iterations, max error = $(@sprintf("%.2e", err))")
            end
            break
        end
        if it == max_iter && settings.verbose
            @warn "Ergodic did not fully converge; max error = $err"
        end
    end
    
    dist ./= sum(dist)
    
    # Reshape back to match v (ie, iy, iw, id)
    dist4 = reshape(dist, (nd, nw, ny, ne))
    dist4 = permutedims(dist4, (4, 3, 2, 1))  # → (ne, ny, nw, nd)
    
    # Edge diagnostics
    if settings.verbose
        wL = sum(dist4[:, :, 1, :])
        wR = sum(dist4[:, :, end, :])
        dL = sum(dist4[:, :, :, 1])
        dR = sum(dist4[:, :, :, end])
        println("Edge mass: w_low=$(round(wL,digits=4)), w_high=$(round(wR,digits=4)), d_low=$(round(dL,digits=4)), d_high=$(round(dR,digits=4))")
    end
    
    return dist4::Array{Float64,4}
end
