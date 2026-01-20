function valfun_noadjust(pea::Vector{Float64}; grid_builder::Function=makegrids)
    grids = grid_builder(pea)
    util, idp_map, d_next_vec = utility_noadjust(grids, pea)

    v    = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    vnew = similar(v)

    gidx = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)
    )
    old = dtp.Ipol(
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(Int, sz.ne, sz.ny, sz.nw, sz.nd)
    )

    pol = dtp.Pol(
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(sz.ne, sz.ny, sz.nw, sz.nd),
        zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    )

    ty = grids.ty

    te_colstoch = is_col_stoch(grids.te)
    te_rowstoch = is_row_stoch(grids.te)
    if !(te_colstoch || te_rowstoch)
        error("grids.te is neither row- nor column-stochastic (numerically).")
    end

    policy_streak = 0

    for iter in 1:sz.maxiter
        queue = similar(v)
        @inbounds for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, iep in 1:sz.ne
            acc = 0.0
            for iyp in 1:sz.ny
                acc += ty[iy, iyp] * v[iep, iyp, iw, id]
            end
            queue[iep, iy, iw, id] = acc
        end

        queuelong = fillin(queue, grids)

        if iter <= sz.earlyiter
            vnew, gidx = maxbellman_noadjust(queuelong, util, idp_map, grids, pea)
        else
            vnew, gidx = tinybellman_noadjust(queuelong, util, idp_map, gidx, grids, pea; te_colstoch=te_colstoch)
        end

        pgap = sum(abs.(gidx.w .- old.w)) + sum(abs.(gidx.s .- old.s))
        policy_streak = (pgap == 0) ? (policy_streak + 1) : 0

        λ = 0.2
        v .= (1 - λ) .* v .+ λ .* vnew

        old.w .= gidx.w; old.d .= gidx.d; old.s .= gidx.s

        if policy_streak >= sz.maxpolit
            pol.w = makepol(gidx.w, grids.wp)
            pol.s = makepol(gidx.s, grids.s)
            @inbounds for id in 1:sz.nd
                pol.d[:, :, :, id] .= d_next_vec[id]
            end
            @inbounds for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
                w_pr = pol.w[ie, iy, iw, id]
                s_pr = pol.s[ie, iy, iw, id]
                E    = grids.ex[ie]
                pol.aa[ie, iy, iw, id] = (1.0 - s_pr) * w_pr
                pol.a[ie, iy, iw, id]  = (s_pr * w_pr) / E
            end
            pol.c = makepol_c(pol.w, pol.d, grids, 0, pea)
            return (v=copy(v), gidx=gidx, pol=pol, g=grids, e=0, idp_map=idp_map, d_next_vec=d_next_vec)
        end
    end

    error("VFI (noadjust) failed")
end
