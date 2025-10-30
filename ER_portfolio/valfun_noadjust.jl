function valfun_noadjust(pea::Vector{Float64})
    beta  = pea[1]
    delta = pea[2]

    # grids first
    grids = makegrids(pea)
    d  = grids.d
    dp = grids.dp
    ddp = (1 - delta) .* d
    iid = [argmin(abs.(dp .- ddp[id])) for id in 1:sz.nd]
    ut  = utility_noadjust(grids, pea)
    tmat = grids.t

    v    = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    vnew = similar(v)

    gidx = dtp.Ipol(
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)),  # a′  (foreign index)
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)),  # aa′ (local  index)
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd))   # d′
    )
    pol = dtp.Pol(
        zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd),        # a′  (foreign level)
        zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd),        # aa′ (local   level)
        zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd),        # d′
        zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)         # c
    )
    old = dtp.Ipol(
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)),
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)),
        Int.(zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd))
    )

    errcode = 0
    in_a_row     = zeros(sz.maxpolit, 1)
    sum_in_a_row = 1000
    gap          = 1000.0
    pgap         = 1000

    Threads.@threads for iter in 1:sz.maxiter
        # queue has same size as v
        queue = zeros(size(v))

        # reshape to multiply only along the first dim (z ≡ (e,y))
        reshaped_v = reshape(v, sz.ne*sz.ny, sz.na, sz.na, sz.nd)
        qq = similar(reshaped_v)

        for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na
            qq[:, iaa, ia, id] = tmat * reshaped_v[:, iaa, ia, id]
        end

        queue = reshape(qq, sz.ne, sz.ny, sz.na, sz.na, sz.nd)

        # interpolate to (npa,npa,npd) if needed
        if sz.na == sz.npa && sz.nd == sz.npd
            queuelong = queue
        else
            queuelong = fillin(queue, grids)  
        end

        # maximize
        if iter <= sz.earlyiter || sum_in_a_row == 0
            if sum_in_a_row > 0
                vnew, gidx = maxbellman_noadjust(queuelong, ut, iid, beta)   
           # else
            #    vnew, gidx = howard_noadjust(queuelong, ut, iid, gidx, beta) 
            end
        else
            vnew, gidx = tinybellman_noadjust(queuelong, ut, iid, gidx, beta) 
        end

        # enforce d′ = iid[id] on both asset dims
        @Threads.threads for id in 1:sz.nd
            @Threads.threads for ia in 1:sz.na
                @Threads.threads for iaa in 1:sz.na
                    @Threads.threads for iy in 1:sz.ny
                        @Threads.threads for ie in 1:sz.ne
                            gidx.d[ie,iy,iaa,ia,id] = iid[id]
                        end
                    end
                end
            end
        end

        # update
        vdiff = vnew - v
        pgap = sum(abs.(gidx.a  - old.a)) +
               sum(abs.(gidx.aa - old.aa)) +
               sum(abs.(gidx.d  - old.d))
        adjfactor = 0.5 * ((beta / (1.0 - beta)) * minimum(vdiff) +
                           (beta / (1.0 - beta)) * maximum(vdiff))
        gap = maximum(abs.(vdiff))
        in_a_row[2:end] = in_a_row[1:end-1]
        in_a_row[1] = pgap
        sum_in_a_row = sum(in_a_row)

        v = vnew .+ adjfactor
        old.a  = gidx.a
        old.aa = gidx.aa
        old.d  = gidx.d

        if mod(iter, 50) == 0 && settings.verbose
            println("iter = ", iter, " function diff = ", gap, " policy diff = ", pgap)
        end
        if iter == sz.maxiter
            errcode = 5
        elseif gap > 10e12
            errcode = 2
        end

        if gap < sz.toler || errcode > 0
            if errcode > 0
                vnew .= 0.0
                gidx.a  .= 0
                gidx.aa .= 0
                gidx.d  .= 0
                pol.a  .= 0.0
                pol.aa .= 0.0
                pol.d  .= 0.0
                pol.c  .= 0.0
                break
            else
                pol.a  = makepol(gidx.a,  grids.ap)    # foreign
                pol.aa = makepol(gidx.aa, grids.aap)   # local
                pol.d  = makepol_d_na(grids.d, delta) 
                pol.c  = makepol_c_twoasset(pol.aa, pol.a, pol.d, grids, 0, pea)  # (aa,a,d), non-adjust
                break
            end
        end
    end

    return (v=vnew::Array{Float64}, gidx=gidx, pol=pol, g=grids::NamedTuple, e=errcode::Int)
end
