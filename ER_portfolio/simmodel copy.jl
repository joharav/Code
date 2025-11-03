function simmodel(answ::NamedTuple)
    v           = answ.v
    pol         = answ.pol            # expects a, aa, d, c
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d

    exg = grids.ex
    yg  = grids.y

    allv        = zeros(sz.nYears, sz.nFirms)
    alla        = zeros(sz.nYears, sz.nFirms)     # foreign
    allaa       = zeros(sz.nYears, sz.nFirms)     # local
    alle        = zeros(sz.nYears, sz.nFirms)
    ally        = zeros(sz.nYears, sz.nFirms)
    alld        = zeros(sz.nYears, sz.nFirms)
    alld_adjust = zeros(sz.nYears, sz.nFirms)
    allc        = zeros(sz.nYears, sz.nFirms)

    # joint (e,y) draws
    phatcdf = cumsum(tmat, dims=2)
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    ls = zeros(sz.nYears + 1, sz.nFirms)
    for ifi in 1:sz.nFirms
        gap = globals.draws[1, ifi] .- cdf_wgt
        ls[1, ifi] = findfirst(gap .< 0.0)
    end

    # initial indices
    astart  = globals.draws[1, :]
    aastart = globals.draws[3, :]  # NEW: local asset init
    dstart  = globals.draws[2, :]

    Threads.@threads for ifi in 1:sz.nFirms
        picka  = min(Int(floor(sz.na * astart[ifi]))  + 1, sz.na)
        pickaa = min(Int(floor(sz.na * aastart[ifi])) + 1, sz.na)
        pickd  = min(Int(floor(sz.nd * dstart[ifi]))  + 1, sz.nd)

        pickey = Int(ls[1, ifi])
        picke = div(pickey - 1, sz.ny) + 1
        picky = mod(pickey - 1, sz.ny) + 1

        vold = v[picke, picky, pickaa, picka, pickd]
        aold  = pol.a[picke, picky, pickaa, picka, pickd]
        aaold = pol.aa[picke, picky, pickaa, picka, pickd]
        dold  = pol.d[picke, picky, pickaa, picka, pickd]
        cold  = pol.c[picke, picky, pickaa, picka, pickd]
        d_adjustold = d_adjust[picke, picky, pickaa, picka, pickd]

        for iti in 1:sz.nYears
            eold = exg[picke]
            yold = yg[picky]

            vprime        = interpol(eold, yold, aaold, aold, dold, grids, v)
            aprime        = interpol(eold, yold, aaold, aold, dold, grids, pol.a)
            aaprime       = interpol(eold, yold, aaold, aold, dold, grids, pol.aa)
            dprime        = interpol(eold, yold, aaold, aold, dold, grids, pol.d)
            d_adjustprime = interpol(eold, yold, aaold, aold, dold, grids, d_adjust)
            cprime        = interpol(eold, yold, aaold, aold, dold, grids, pol.c)

            gap = globals.draws[iti + 1, ifi] .- phatcdf[Int(ls[iti, ifi]), :]
            ls[iti + 1, ifi] = findfirst(gap .< 0.0)

            pickey_next = Int(ls[iti + 1, ifi])
            picke = div(pickey_next - 1, sz.ny) + 1
            picky = mod(pickey_next - 1, sz.ny) + 1

            eprime = exg[picke]
            yprime = yg[picky]

            allv[iti, ifi]        = vprime[1]
            alla[iti, ifi]        = aprime[1]
            allaa[iti, ifi]       = aaprime[1]
            alle[iti, ifi]        = eprime
            ally[iti, ifi]        = yprime
            alld[iti, ifi]        = dprime[1]
            alld_adjust[iti, ifi] = d_adjustprime[1]
            allc[iti, ifi]        = cprime[1]

            # update states
            vold  = vprime[1]
            aold  = aprime[1]
            aaold = aaprime[1]
            dold  = dprime[1]
            d_adjustold = d_adjustprime[1]
            cold  = cprime[1]
        end
    end

    adjust_indicator = zeros(size(alld))
    for i in 1:sz.nYears, j in 1:sz.nFirms
        adjust_indicator[i, j] = (alld_adjust[i, j] == alld[i, j]) ? 1.0 : 0.0
    end

    return (
        v = allv,
        d = alld,
        a = alla,
        aa = allaa,
        ex = alle,
        y = ally,
        d_adjust = alld_adjust,
        adjust_indicator = adjust_indicator,
        c = allc
    )
end
