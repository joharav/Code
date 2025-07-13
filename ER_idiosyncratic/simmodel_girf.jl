function simmodel_girf(answ::NamedTuple, T_shock::Int)
    v           = answ.v
    pol         = answ.pol
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d

    exg = grids.ex
    yg  = grids.y

    nT, nI = sz.nYears, sz.nFirms
    allv, alla, alle, ally = zeros(nT, nI), zeros(nT, nI), zeros(nT, nI), zeros(nT, nI)
    alld, alld_adjust, allc = zeros(nT, nI), zeros(nT, nI), zeros(nT, nI)

    phatcdf = cumsum(tmat, dims=2)
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])  # assumes ex & y jointly represented
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

    astart, dstart = globals.draws[1, :], globals.draws[2, :]
    shock_val = maximum(exg)

    Threads.@threads for i in 1:nI
        picka = min(Int(floor(sz.na * astart[i])) + 1, sz.na)
        pickd = min(Int(floor(sz.nd * dstart[i])) + 1, sz.nd)

        picke = all_picke[1, i]
        picky = all_picky[1, i]

        vold = v[picke, picky, picka, pickd]
        aold = pol.a[picke, picky, picka, pickd]
        dold = pol.d[picke, picky, picka, pickd]
        cold = pol.c[picke, picky, picka, pickd]
        d_adjustold = d_adjust[picke, picky, picka, pickd]

        for t in 1:nT
            eold = (t == T_shock) ? shock_val : exg[picke]
            yold = yg[picky]

            aprime = interpol(eold, yold, aold, dold, grids, pol.a)
            dprime = interpol(eold, yold, aold, dold, grids, pol.d)
            d_adjustprime = interpol(eold, yold, aold, dold, grids, d_adjust)
            cprime = interpol(eold, yold, aold, dold, grids, pol.c)
            vprime = interpol(eold, yold, aold, dold, grids, v)

            picke = all_picke[t + 1, i]
            picky = all_picky[t + 1, i]

            alle[t, i] = (t == T_shock) ? shock_val : exg[picke]
            ally[t, i] = yg[picky]
            allv[t, i] = vprime[1]
            alla[t, i] = aprime[1]
            alld[t, i] = dprime[1]
            alld_adjust[t, i] = d_adjustprime[1]
            allc[t, i] = cprime[1]

            # update states
            vold, aold, dold = vprime[1], aprime[1], dprime[1]
            d_adjustold, cold = d_adjustprime[1], cprime[1]
        end
    end

    adjust_indicator = map((x, y) -> x == y ? 1.0 : 0.0, alld_adjust, alld)

    return (
        v = allv, d = alld, a = alla, ex = alle, y = ally,
        d_adjust = alld_adjust, adjust_indicator = adjust_indicator, c = allc
    )
end
