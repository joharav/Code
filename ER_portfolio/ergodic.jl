function compute_ergodic(ans::NamedTuple; tol=sz.distol, max_iter=sz.maxditer)
    g    = ans.g
    tmat = g.t  # size (ne*ny)×(ne*ny)

    ne, ny, na, nd = sz.ne, sz.ny, sz.na, sz.nd
    nzy     = ne * ny
    nstates = ne * ny * na * na * nd

    # Policies must be 5-D: (ie, iy, iaa, ia, id)
    pol_a  = ans.pol.a
    pol_aa = ans.pol.aa
    pol_d  = ans.pol.d

    @assert size(pol_a)  == (ne, ny, na, na, nd)
    @assert size(pol_aa) == (ne, ny, na, na, nd)
    @assert size(pol_d)  == (ne, ny, na, na, nd)
    @assert size(tmat)   == (nzy, nzy)

    # Linear index consistent with reshape(..., (nd,na,na,ny,ne)) in column-major Julia
    flat_idx(ie, iy, iaa, ia, id) =
        (((((ie-1)*ny + (iy-1))*na + (iaa-1))*na + (ia-1))*nd + (id-1)) + 1

    find_index(val, grid) = clamp(searchsortedfirst(grid, val), 1, length(grid))

    dist     = fill(1.0 / nstates, nstates)
    dist_new = similar(dist)

    for it in 1:max_iter
        fill!(dist_new, 0.0)

        @inbounds for ie in 1:ne, iy in 1:ny, iaa in 1:na, ia in 1:na, id in 1:nd
            idx  = flat_idx(ie, iy, iaa, ia, id)
            mass = dist[idx]
            mass == 0.0 && continue

            # Next-period policy indices
            ia′  = find_index(pol_a[ie, iy, iaa, ia, id], g.a)
            iaa′ = find_index(pol_aa[ie, iy, iaa, ia, id], g.aa)
            id′  = find_index(pol_d[ie, iy, iaa, ia, id], g.d)

            # Transition over z = (e,y)
            iz = (iy - 1) * ne + ie
            for iz′ in 1:nzy
                ie′ = mod1(iz′, ne)
                iy′ = div(iz′ - 1, ne) + 1
                prob = tmat[iz′, iz]
                prob == 0.0 && continue
                idx′ = flat_idx(ie′, iy′, iaa′, ia′, id′)
                dist_new[idx′] += prob * mass
            end
        end

        err = maximum(abs.(dist_new .- dist))
        dist, dist_new = dist_new, dist
        if err < tol
            # println("Ergodic distribution converged after $it iterations, max error = $err")
            break
        end
        if it == max_iter && settings.verbose
            @warn "Ergodic distribution did not fully converge; max error = $err"
        end
    end

    dist ./= sum(dist)

    # Reshape back to match v (ie, iy, iaa, ia, id)
    dist5 = reshape(dist, (nd, na, na, ny, ne))
    dist5 = permutedims(dist5, (5, 4, 3, 2, 1))  # → (ne, ny, na, na, nd)

    # Quick edge diagnostics (optional)
    if settings.verbose
        aL  = sum(dist5[:,:,:,1,:]);  aR  = sum(dist5[:,:,: ,end,:])
        aaL = sum(dist5[:,:,1,:,:]);  aaR = sum(dist5[:,: ,end,:,:])
        dL  = sum(dist5[:,:,:,:,1]);  dR  = sum(dist5[:,:,:,:,end])
        tot_edge = aL+aR+aaL+aaR+dL+dR
        println("Edge mass (a,aa,d) = $(round(tot_edge,digits=4))")
    end

    return dist5::Array{Float64,5}
end
