function compute_ergodic(ans::NamedTuple; tol=sz.distol, max_iter=sz.maxditer)
    g    = ans.g
    tmat = g.t                      # (ne*ny)×(ne*ny)
    na   = sz.na; nd = sz.nd; ne = sz.ne; ny = sz.ny
    nzy  = ne * ny
    nstates = na * na * nd * nzy

    pol_a  = ans.pol.a   # [ne,ny,iaa,ia,id] levels
    pol_aa = ans.pol.aa  # [ne,ny,iaa,ia,id] levels
    pol_d  = ans.pol.d   # [ne,ny,iaa,ia,id] levels

    # flatten (ie,iy,id,iaa,ia)
    flat_idx(ie,iy,id,iaa,ia) = (((((ie-1)*ny + (iy-1))*nd + (id-1))*na + (iaa-1))*na + ia) + 1

    find_index(val, grid) = clamp(searchsortedfirst(grid, val), 1, length(grid))

    dist     = fill(1.0 / nstates, nstates)
    dist_new = similar(dist)

    for it in 1:max_iter
        fill!(dist_new, 0.0)
        for ie in 1:ne, iy in 1:ny, id in 1:nd, iaa in 1:na, ia in 1:na
            idx  = flat_idx(ie,iy,id,iaa,ia)
            mass = dist[idx]

            aa1 = pol_aa[ie,iy,iaa,ia,id]
            a1  = pol_a[ie,iy,iaa,ia,id]
            d1  = pol_d[ie,iy,iaa,ia,id]

            iaa′ = find_index(aa1, g.aa)
            ia′  = find_index(a1,  g.a)
            id′  = find_index(d1,  g.d)

            iz   = (iy - 1) * ne + ie
            for iz′ in 1:nzy
                ie′ = mod1(iz′, ne)
                iy′ = div(iz′ - 1, ne) + 1
                prob = tmat[iz′, iz]
                idx′ = flat_idx(ie′,iy′,id′,iaa′,ia′)
                dist_new[idx′] += prob * mass
            end
        end
        err = maximum(abs.(dist_new .- dist))
        if err < tol
            println("Ergodic distribution converged after $it iterations, max error = $err")
            break
        end
        dist, dist_new = dist_new, dist
    end

    dist ./= sum(dist)
    dist5d = reshape(dist, (na, na, nd, ny, ne))
    dist5d = permutedims(dist5d, (5, 4, 3, 1, 2))  # [ie, iy, id, iaa, ia]

    # edge mass checks
    aL_left   = sum(dist5d[:,:,: , 1, :])        # aa_min
    aL_right  = sum(dist5d[:,:,: , end, :])      # aa_max
    aF_left   = sum(dist5d[:,:,: , :, 1])        # a_min
    aF_right  = sum(dist5d[:,:,: , :, end])      # a_max
    d_left    = sum(dist5d[:,:, 1, :, :])
    d_right   = sum(dist5d[:,:, end, :, :])

    if (aL_left + aL_right + aF_left + aF_right + d_left + d_right) > 0.001
        println("%%%%%%%% WARNING %%%%%%%%%%")
        println("High ergodic mass on grid edges (two-asset). Consider expanding grids.")
        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    else
        println("→ Grid coverage looks adequate.")
        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    end

    return dist5d::Array{Float64,5}
end
