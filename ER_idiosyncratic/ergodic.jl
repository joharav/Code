function compute_ergodic(ans::NamedTuple; tol=sz.distol, max_iter=sz.maxditer)
    g = ans.g
    tmat = g.t  # full transition over z = (ie, iy) flattened: size (nzy, nzy)

    na, nd, ne, ny = sz.na, sz.nd, sz.ne, sz.ny
    nzy = ne * ny
    nstates = na * nd * ne * ny

    pol_a = ans.pol.a  # size [ne, ny, nd, na]
    pol_d = ans.pol.d  # size [ne, ny, nd, na]

    # Flatten index
    function flat_idx(ie, iy, id, ia)
        return (((ie-1) * ny + (iy-1)) * nd + (id-1)) * na + ia
    end

    function find_index(val, grid)
        return clamp(searchsortedfirst(grid, val), 1, length(grid))
    end

    dist = fill(1.0 / nstates, nstates)
    dist_new = similar(dist)

    for it in 1:max_iter
        fill!(dist_new, 0.0)

        for ie in 1:ne, iy in 1:ny, id in 1:nd, ia in 1:na
            idx = flat_idx(ie, iy, id, ia)
            mass = dist[idx]

            a_next = pol_a[ie, iy, id, ia]
            d_next = pol_d[ie, iy, id, ia]

            ia′ = find_index(a_next, g.a)
            id′ = find_index(d_next, g.d)

            # current flat z index and iterate over all future states
            iz = (iy - 1) * ne + ie
            for iz′ in 1:nzy
                ie′ = mod1(iz′, ne)
                iy′ = div(iz′ - 1, ne) + 1

                prob = tmat[iz′, iz]
                idx′ = flat_idx(ie′, iy′, id′, ia′)
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

    dist = dist ./ sum(dist)

    dist4d = reshape(dist, (na, nd, ny, ne))
    dist4d = permutedims(dist4d, (4, 3, 2, 1))  # [ie, iy, id, ia]

    # === Check for endpoint mass on assets ===
    a_left_wt  = sum(dist4d[:, :, :, 1])      # a_min
    a_right_wt = sum(dist4d[:, :, :, end])    # a_max

    # === Check for endpoint mass on durables ===
    d_left_wt  = sum(dist4d[:, :, 1, :])      # d_min
    d_right_wt = sum(dist4d[:, :, end, :])    # d_max

    total_edge_wt = a_left_wt + a_right_wt + d_left_wt + d_right_wt
    total_edge_a = a_left_wt + a_right_wt 
    total_edge_d = d_left_wt + d_right_wt

    if total_edge_wt > 0.001 ||  total_edge_a > 0.001 || total_edge_d > 0.001# or your preferred threshold
        println("%%%%%%%% WARNING %%%%%%%%%%")
        println("High ergodic mass on grid edges:")
        println("Assets → Left = $(round(a_left_wt, digits=4)), Right = $(round(a_right_wt, digits=4)), Total = $(round(total_edge_a, digits=4))")
        println("Durables → Left = $(round(d_left_wt, digits=4)), Right = $(round(d_right_wt, digits=4)), Total = $(round(total_edge_d, digits=4))")
        println("→ Suggest expanding grids (assets or durables).")
        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%")

       
    else
        println("→ Grid coverage looks adequate.")
        println("%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    end


    return dist4d::Array{Float64, 4}

    
end
