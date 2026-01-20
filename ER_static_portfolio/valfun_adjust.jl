# ==========================================================================
# 4D MODEL: Value Function Iteration for ADJUSTMENT regime
# State: (e, y, w, d)
# Policy: (w', d', s*)
# ==========================================================================

function valfun_adjust(pea::Vector{Float64})
    beta = pea[1]

    # Build grids
    grids = grid_builder(pea)
    ut = utility(grids, pea)
    tmat_joint = grids.t      # joint (e,y) transition
    
    # Initialize value function and policies
    v = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    vnew = similar(v)

    gidx = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),  # w'
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),  # d'
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)   # s
    )
    old = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)
    )

    pol = dtp.Pol(
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # w
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # d
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # s
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # c
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),   # aa (derived)
        zeros(sz.ne, sz.ny, sz.nw, sz.nd)    # a (derived)
    )

    errcode = 0
    gap = Inf
    pgap = typemax(Int)
    policy_streak = 0

    for iter in 1:sz.maxiter
        # 1) Compute expected continuation value E[V(e',y',w,d) | e,y]
        # This integrates over (e',y') transitions
        # Note: In 4D model, the w' -> w_{t+1} mapping depends on portfolio choice,
        # so we can't fully pre-compute E[V]. Instead, we pass V to maxbellman
        # which handles the portfolio-dependent wealth evolution.
        
        # For the joint (e,y) transition, reshape and multiply
        reshaped_v = reshape(v, sz.ne * sz.ny, sz.nw, sz.nd)
        qq = similar(reshaped_v)
        
        @inbounds for id in 1:sz.nd, iw in 1:sz.nw
            # E[V(e',y',w,d) | e,y] for fixed (w,d)
            qq[:, iw, id] = tmat_joint * reshaped_v[:, iw, id]
        end
        
        queue = reshape(qq, sz.ne, sz.ny, sz.nw, sz.nd)

        # 2) Bellman update with integrated portfolio choice
        if iter <= sz.earlyiter
            vnew, gidx = maxbellman(queue, ut, grids, pea)
        else
            vnew, gidx = tinybellman(queue, ut, gidx, grids, pea)
        end

        # 3) Convergence diagnostics
        vdiff = vnew .- v
        gap = maximum(abs.(vdiff))

        pgap = sum(abs.(gidx.w .- old.w)) +
               sum(abs.(gidx.d .- old.d)) +
               sum(abs.(gidx.s .- old.s))

        if pgap == 0
            policy_streak += 1
        else
            policy_streak = 0
        end

        # 4) Update value with acceleration
        adjfactor = 0.5 * ((beta / (1.0 - beta)) * minimum(vdiff) +
                          (beta / (1.0 - beta)) * maximum(vdiff))
        v = vnew .+ adjfactor

        # Store old policies
        old.w .= gidx.w
        old.d .= gidx.d
        old.s .= gidx.s

        if settings.verbose && (iter % 50 == 0)
            println("Adjust iter=$iter gap=$(@sprintf("%.2e", gap)) pgap=$pgap streak=$policy_streak")
        end

        # Failure checks
        if gap > 1e13
            errcode = 2
            break
        elseif iter == sz.maxiter
            errcode = 5
            break
        end

        # Success: stable policy
        if policy_streak >= sz.maxpolit
            # Extract policy levels from indices
            pol.w = makepol(gidx.w, grids.wp)
            pol.d = makepol(gidx.d, grids.dp)
            pol.s = makepol(gidx.s, grids.s)
            
            # Derive aa and a from w and s
            # aa = (1-s) * w', a = s * w' / E
            for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
                w_pr = pol.w[ie, iy, iw, id]
                s_pr = pol.s[ie, iy, iw, id]
                E = grids.ex[ie]
                
                pol.aa[ie, iy, iw, id] = (1.0 - s_pr) * w_pr
                pol.a[ie, iy, iw, id] = (s_pr * w_pr) / E
            end
            
            # Compute implied consumption
            pol.c = makepol_c(pol.w, pol.d, pol.s, grids, 1, pea)
            break
        end
    end

    if errcode > 0
        error("VFI (adjust) failed: errcode=$errcode gap=$gap pgap=$pgap")
    end

    return (v=vnew, gidx=gidx, pol=pol, g=grids, e=errcode)
end


# ==========================================================================
# Local search (tinybellman) for adjustment regime
# ==========================================================================

function tinybellman(queue::Array{Float64,4}, util::Array{Float64,6},
                    old_gidx::dtp.Ipol, grids::NamedTuple, pea::Vector{Float64})
    
    beta = pea[1]
    rr = (1 / beta) - 1
    rr_star = pea[9]
    kappa = pea[11]
    
    e_grid = grids.ex
    w_grid = grids.w
    wp_grid = grids.wp
    s_grid = grids.s
    trans_e = grids.te
    
    vnew = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    gidx = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)
    )

    Threads.@threads for id in 1:sz.nd
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    # Local search around previous optimum
                    old_iwp = old_gidx.w[ie, iy, iw, id]
                    old_idp = old_gidx.d[ie, iy, iw, id]
                    old_is = old_gidx.s[ie, iy, iw, id]
                    
                    # Search ranges
                    iwp_lo = max(1, old_iwp - sz.pad)
                    iwp_hi = min(sz.npw, old_iwp + sz.pad)
                    idp_lo = max(1, old_idp - sz.pad)
                    idp_hi = min(sz.npd, old_idp + sz.pad)
                    is_lo = max(1, old_is - 3)
                    is_hi = min(sz.ns, old_is + 3)
                    
                    vstar = -Inf
                    wstar, dstar, sstar = old_iwp, old_idp, old_is
                    E_now = e_grid[ie]
                    
                    for idp in idp_lo:idp_hi
                        for iwp in iwp_lo:iwp_hi
                            u_flow = util[ie, iy, iw, id, iwp, idp]
                            
                            if u_flow > -1e9
                                w_next = wp_grid[iwp]
                                
                                # Local search for s
                                for is in is_lo:is_hi
                                    s = s_grid[is]
                                    trans_cost = kappa * s * w_next
                                    
                                    EV = 0.0
                                    for ie_next in 1:sz.ne
                                        E_next = e_grid[ie_next]
                                        w_realized = (1.0 - s) * w_next * (1.0 + rr) + 
                                                    s * w_next * (1.0 + rr_star) * (E_next / E_now) -
                                                    trans_cost
                                        w_realized = clamp(w_realized, w_grid[1], w_grid[end])
                                        
                                        iw_L, iw_U, wt = bracket_grid(w_realized, w_grid)
                                        V_interp = (1.0 - wt) * queue[ie_next, iy, iw_L, idp] +
                                                   wt * queue[ie_next, iy, iw_U, idp]
                                        EV += trans_e[ie, ie_next] * V_interp
                                    end
                                    
                                    bellman = u_flow + beta * EV
                                    if bellman > vstar
                                        vstar = bellman
                                        wstar = iwp
                                        dstar = idp
                                        sstar = is
                                    end
                                end
                            end
                        end
                    end
                    
                    vnew[ie, iy, iw, id] = vstar
                    gidx.w[ie, iy, iw, id] = wstar
                    gidx.d[ie, iy, iw, id] = dstar
                    gidx.s[ie, iy, iw, id] = sstar
                end
            end
        end
    end
    
    return vnew, gidx
end
