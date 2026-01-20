# ==========================================================================
# Stationary distribution for 4D model
# State: (e, y, w, d)
# Uses:
#   - grids.t : joint transition over combined (e,y) index of length ne*ny
#   - answ.pol.w, answ.pol.d, answ.pol.s : policy levels on state grid
# Wealth transition uses portfolio returns + FX, then interpolates onto w-grid
# Durable transition interpolates onto d-grid (or collapses if exact)
# ==========================================================================

@inline function ey_index(ie::Int, iy::Int, ny::Int)
    return (ie - 1) * ny + iy
end

@inline function ey_decompose(j::Int, ny::Int)
    ie = div(j - 1, ny) + 1
    iy = mod(j - 1, ny) + 1
    return ie, iy
end

@inline function bracket_grid(x::Float64, g::AbstractVector{<:Real})
    n = length(g)
    if x <= g[1]
        return 1, 1, 0.0
    elseif x >= g[n]
        return n, n, 0.0
    else
        j = searchsortedlast(g, x)
        xL = g[j]
        xU = g[j+1]
        wt = (xU == xL) ? 0.0 : (x - xL) / (xU - xL)
        return j, j+1, wt
    end
end

function stationary_dist_4d(answ::NamedTuple;
                            maxiter::Int = 5_000,
                            tol::Float64 = 1e-10)

    grids = answ.g
    tmat  = grids.t                 # (ne*ny) x (ne*ny)
    exg   = grids.ex
    yg    = grids.y
    wg    = grids.w
    dg    = grids.d
    sg    = grids.s

    # parameters for wealth realization
    beta    = answ.pea[1]
    rr      = (1 / beta) - 1
    rr_star = answ.pea[9]
    kappa   = answ.pea[11]

    ne, ny, nw, nd = sz.ne, sz.ny, sz.nw, sz.nd

    # policies on state grid (levels)
    pol_w = answ.pol.w
    pol_d = answ.pol.d
    pol_s = answ.pol.s

    # initialize uniform distribution
    μ  = fill(1.0 / (ne * ny * nw * nd), ne, ny, nw, nd)
    μn = similar(μ)

    # precompute endogenous mapping weights for each (ie,iy,iw,id) to w-grid and d-grid
    wL = Array{Int}(undef, ne, ny, nw, nd)
    wU = Array{Int}(undef, ne, ny, nw, nd)
    ww = Array{Float64}(undef, ne, ny, nw, nd)

    dL = Array{Int}(undef, ne, ny, nw, nd)
    dU = Array{Int}(undef, ne, ny, nw, nd)
    dw = Array{Float64}(undef, ne, ny, nw, nd)

    @inbounds for id in 1:nd, iw in 1:nw, iy in 1:ny, ie in 1:ne
        # chosen next-period (end of period) targets
        w_next = pol_w[ie, iy, iw, id]
        d_next = pol_d[ie, iy, iw, id]

        iL, iU, wt = bracket_grid(w_next, wg)
        wL[ie,iy,iw,id] = iL
        wU[ie,iy,iw,id] = iU
        ww[ie,iy,iw,id] = wt

        jL, jU, wt2 = bracket_grid(d_next, dg)
        dL[ie,iy,iw,id] = jL
        dU[ie,iy,iw,id] = jU
        dw[ie,iy,iw,id] = wt2
    end

    # iterate μ' = P' μ
    for it in 1:maxiter
        fill!(μn, 0.0)

        @inbounds for id in 1:nd, iw in 1:nw, iy in 1:ny, ie in 1:ne
            mass = μ[ie,iy,iw,id]
            mass == 0.0 && continue

            # exogenous transition over (e,y)
            row = ey_index(ie, iy, ny)

            # portfolio-induced wealth realization depends on next e'
            # we will handle e' inside the ey-loop, then interpolate wealth back to w-grid.

            # fixed chosen portfolio share s at this state
            s = pol_s[ie,iy,iw,id]
            s = clamp(s, 0.0, 1.0)

            # chosen end-of-period saved wealth (before returns)
            w_next = pol_w[ie,iy,iw,id]

            # durable interpolation weights (do not depend on e')
            jL = dL[ie,iy,iw,id]
            jU = dU[ie,iy,iw,id]
            wt_d = dw[ie,iy,iw,id]

            E_now = exg[ie]
            trans_cost = kappa * s * w_next

            for col in 1:(ne*ny)
                p = tmat[row, col]
                p == 0.0 && continue

                iep, iyp = ey_decompose(col, ny)
                E_next = exg[iep]

                # realized wealth in pesos next period
                w_real = (1.0 - s) * w_next * (1.0 + rr) +
                         s * w_next * (1.0 + rr_star) * (E_next / max(E_now, 1e-12)) -
                         trans_cost

                w_real = min(max(w_real, wg[1]), wg[end])
                iL, iU, wt_w = bracket_grid(w_real, wg)

                pmass = mass * p

                # bilinear distribution across (w,d)
                # weights:
                #   w: (1-wt_w), wt_w
                #   d: (1-wt_d), wt_d
                μn[iep, iyp, iL, jL] += pmass * (1-wt_w) * (1-wt_d)
                μn[iep, iyp, iU, jL] += pmass * (wt_w)   * (1-wt_d)
                μn[iep, iyp, iL, jU] += pmass * (1-wt_w) * (wt_d)
                μn[iep, iyp, iU, jU] += pmass * (wt_w)   * (wt_d)
            end
        end

        # normalize (numerical drift)
        sμ = sum(μn)
        μn ./= sμ

        err = maximum(abs.(μn .- μ))
        μ, μn = μn, μ
        if err < tol
            break
        end
    end

    return μ
end
