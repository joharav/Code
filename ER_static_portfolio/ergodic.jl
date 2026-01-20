# ==========================================================================
# 4D MODEL: Compute ergodic distribution
# State: (e, y, w, d)
# Uses bilinear mass-splitting on (w,d) given policy levels pol_w, pol_d
# ==========================================================================

function compute_ergodic(ans::NamedTuple; tol=sz.distol, max_iter=sz.maxditer)
    g = ans.g
    tmat = g.t  # (ne*ny) × (ne*ny)

    ne, ny, nw, nd = sz.ne, sz.ny, sz.nw, sz.nd
    nzy = ne * ny
    nstates = ne * ny * nw * nd

    pol_w = ans.pol.w   # levels (ne,ny,nw,nd)
    pol_d = ans.pol.d   # levels (ne,ny,nw,nd)

    @assert size(pol_w) == (ne, ny, nw, nd)
    @assert size(pol_d) == (ne, ny, nw, nd)
    @assert size(tmat) == (nzy, nzy)

    # Linear index (column-major; fastest-changing index is id)
    @inline flat_idx(ie, iy, iw, id) =
        ((((ie-1)*ny + (iy-1))*nw + (iw-1))*nd + (id-1)) + 1

    # Initialize uniform distribution on full state space
    dist = fill(1.0 / nstates, nstates)
    dist_new = similar(dist)

    # Precompute z-index mapping (consistent with your tmat usage)
    @inline z_idx(ie, iy) = (iy - 1) * ne + ie  # 1..ne*ny
    @inline ie_from_z(iz) = mod1(iz, ne)
    @inline iy_from_z(iz) = div(iz - 1, ne) + 1

    for it in 1:max_iter
        fill!(dist_new, 0.0)

        @inbounds for ie in 1:ne, iy in 1:ny, iw in 1:nw, id in 1:nd
            idx = flat_idx(ie, iy, iw, id)
            mass = dist[idx]
            mass == 0.0 && continue

            # Policy-implied next-period levels
            w′ = pol_w[ie, iy, iw, id]
            d′ = pol_d[ie, iy, iw, id]

            # Bracket on STATE grids (value-based; correct even if policies came from dp/wp grids)
            iwL, iwU, ww = brack1d(g.w, w′)
            idL, idU, wd = brack1d(g.d, d′)

            # Weights for 2D bilinear split
            w00 = (1.0 - ww) * (1.0 - wd)
            w10 = (      ww) * (1.0 - wd)
            w01 = (1.0 - ww) * (      wd)
            w11 = (      ww) * (      wd)

            iz = z_idx(ie, iy)

            # Transition over z' with probabilities tmat[z', z]
            for iz′ in 1:nzy
                prob = tmat[iz′, iz]
                prob == 0.0 && continue

                ie′ = ie_from_z(iz′)
                iy′ = iy_from_z(iz′)

                # Push mass to the 4 neighboring (w,d) nodes
                dist_new[flat_idx(ie′, iy′, iwL, idL)] += prob * mass * w00
                dist_new[flat_idx(ie′, iy′, iwU, idL)] += prob * mass * w10
                dist_new[flat_idx(ie′, iy′, iwL, idU)] += prob * mass * w01
                dist_new[flat_idx(ie′, iy′, iwU, idU)] += prob * mass * w11
            end
        end

        # Convergence check
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

    # Normalize (rounding safety)
    s = sum(dist)
    s > 0 || error("compute_ergodic: distribution collapsed to zero mass")
    dist ./= s

    # Reshape back to (ne, ny, nw, nd)
    dist4 = reshape(dist, (nd, nw, ny, ne))
    dist4 = permutedims(dist4, (4, 3, 2, 1))

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
