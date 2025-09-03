function simmodel_girf(answ::NamedTuple, T_shock::Int)
    v           = answ.v
    pol         = answ.pol              # expects pol.a, pol.aa, pol.d, pol.c with dims [ie,iy,iaa,ia,id]
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d

    exg = grids.ex
    yg  = grids.y

    nT, nI = sz.nYears, sz.nFirms
    allv  = zeros(nT, nI)
    alla  = zeros(nT, nI)     # foreign a
    allaa = zeros(nT, nI)     # local  aa
    alle  = zeros(nT, nI)
    ally  = zeros(nT, nI)
    alld  = zeros(nT, nI)
    alld_adjust = zeros(nT, nI)
    allc  = zeros(nT, nI)

    # joint (e,y) draws as before
    phatcdf = cumsum(tmat, dims=2)
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    ls = zeros(nT + 1, nI)
    for i in 1:nI
        gap = globals.draws[1, i] .- cdf_wgt
        ls[1, i] = findfirst(gap .< 0)
    end
    for i in 1:nI, t in 1:nT
        gap = globals.draws[t + 1, i] .- phatcdf[Int(ls[t, i]), :]
        ls[t + 1, i] = findfirst(gap .< 0)
    end
    all_picke = div.(Int.(ls .- 1), sz.ny) .+ 1
    all_picky = mod.(Int.(ls .- 1), sz.ny) .+ 1

    # initial conditions (reuse extra rows for aa)
    astart   = globals.draws[1, :]
    aastart  = globals.draws[3, :]  # NEW: local asset init from row 3
    dstart   = globals.draws[2, :]

    shock_val = maximum(exg)

    Threads.@threads for i in 1:nI
        picka  = min(Int(floor(sz.na * astart[i]))  + 1, sz.na)
        pickaa = min(Int(floor(sz.na * aastart[i])) + 1, sz.na)
        pickd  = min(Int(floor(sz.nd * dstart[i]))  + 1, sz.nd)

        picke = all_picke[1, i]
        picky = all_picky[1, i]

        # levels at t=1
        vold = v[picke, picky, pickaa, picka, pickd]
        aold  = pol.a[picke, picky, pickaa, picka, pickd]
        aaold = pol.aa[picke, picky, pickaa, picka, pickd]
        dold  = pol.d[picke, picky, pickaa, picka, pickd]
        cold  = pol.c[picke, picky, pickaa, picka, pickd]
        d_adjustold = d_adjust[picke, picky, pickaa, picka, pickd]

        for t in 1:nT
            eold = (t == T_shock) ? shock_val : exg[picke]
            yold = yg[picky]

            # 5D interpolation: (e, y, aa, a, d)
            aprime        = interpol(eold, yold, aaold, aold, dold, grids, pol.a)
            aaprime       = interpol(eold, yold, aaold, aold, dold, grids, pol.aa)
            dprime        = interpol(eold, yold, aaold, aold, dold, grids, pol.d)
            d_adjustprime = interpol(eold, yold, aaold, aold, dold, grids, d_adjust)
            cprime        = interpol(eold, yold, aaold, aold, dold, grids, pol.c)
            vprime        = interpol(eold, yold, aaold, aold, dold, grids, v)

            # next (e,y)
            picke = all_picke[t + 1, i]
            picky = all_picky[t + 1, i]

            alle[t, i] = (t == T_shock) ? shock_val : exg[picke]
            ally[t, i] = yg[picky]
            allv[t, i] = vprime[1]
            alla[t, i] = aprime[1]
            allaa[t, i] = aaprime[1]
            alld[t, i] = dprime[1]
            alld_adjust[t, i] = d_adjustprime[1]
            allc[t, i] = cprime[1]

            # update
            vold  = vprime[1]
            aold  = aprime[1]
            aaold = aaprime[1]
            dold  = dprime[1]
            d_adjustold = d_adjustprime[1]
            cold  = cprime[1]
        end
    end

    adjust_indicator = map((x, y) -> x == y ? 1.0 : 0.0, alld_adjust, alld)

    return (
        v = allv, d = alld, a = alla, aa = allaa,
        ex = alle, y = ally, d_adjust = alld_adjust,
        adjust_indicator = adjust_indicator, c = allc
    )
end
