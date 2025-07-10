function simmodel(answ::NamedTuple)
    # Extract inputs
    v           = answ.v
    pol         = answ.pol
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d

    # Grids
    exg = grids.ex
    yg  = grids.y  # NEW

    # Initialize outputs
    allv        = zeros(sz.nYears, sz.nFirms)
    alla        = zeros(sz.nYears, sz.nFirms)
    alle        = zeros(sz.nYears, sz.nFirms)
    ally        = zeros(sz.nYears, sz.nFirms)  # NEW
    alld        = zeros(sz.nYears, sz.nFirms)
    alld_adjust = zeros(sz.nYears, sz.nFirms)
    allc        = zeros(sz.nYears, sz.nFirms)

    # CDF for simulation
    phatcdf = cumsum(tmat, dims=2)
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])  # assumes ex & y jointly represented

    # Shock indices
    ls = zeros(sz.nYears + 1, sz.nFirms)  # location indices
    for ifi in 1:sz.nFirms
        gap = globals.draws[1, ifi] .- cdf_wgt
        gap = gap .< 0.0
        ls[1, ifi] = findfirst(gap)
    end

    # Initial positions
    global astart = globals.draws[1, :]
    global dstart = globals.draws[2, :]
    
    Threads.@threads for ifi in 1:sz.nFirms
        picka = min(Int(floor(sz.na * astart[ifi])) + 1, sz.na)
        pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
        pickey = Int(ls[1, ifi])  # combined index for ex and y

        # Decompose if ex and y are jointly indexed:
        picke = div(pickey - 1, sz.ny) + 1
        picky = mod(pickey - 1, sz.ny) + 1

        vold =    v[picke, picky, picka, pickd]
        aold = pol.a[picke, picky, picka, pickd]
        dold = pol.d[picke, picky, picka, pickd]
        cold = pol.c[picke, picky, picka, pickd]
        d_adjustold = d_adjust[picke, picky, picka, pickd]

        for iti in 1:sz.nYears
            eold = exg[picke]
            yold = yg[picky]

            # Interpolation in 4D: (e, y, a, d)
            vprime        = interpol(eold, yold, aold, dold, grids, v)
            aprime        = interpol(eold, yold, aold, dold, grids, pol.a)
            dprime        = interpol(eold, yold, aold, dold, grids, pol.d)
            d_adjustprime = interpol(eold, yold, aold, dold, grids, d_adjust)
            cprime        = interpol(eold, yold, aold, dold, grids, pol.c)

            # Update shock index
            gap = globals.draws[iti + 1, ifi] .- phatcdf[Int(ls[iti, ifi]), :]
            gap = gap .< 0.0
            ls[iti + 1, ifi] = findfirst(gap)

            # Decompose again
            pickey_next = Int(ls[iti + 1, ifi])
            picke = div(pickey_next - 1, sz.ny) + 1
            picky = mod(pickey_next - 1, sz.ny) + 1

            eprime = exg[picke]
            yprime = yg[picky]

            # Store
            allv[iti, ifi]           = vprime[1]
            alla[iti, ifi]           = aprime[1]
            alle[iti, ifi]           = eprime
            ally[iti, ifi]           = yprime
            alld[iti, ifi]           = dprime[1]
            alld_adjust[iti, ifi]    = d_adjustprime[1]
            allc[iti, ifi]           = cprime[1]

            # Update states
            vold = vprime[1]
            aold = aprime[1]
            dold = dprime[1]
            d_adjustold = d_adjustprime[1]
            cold = cprime[1]
        end
    end

    # Adjustment indicator
    adjust_indicator = zeros(size(alld))
    for i in 1:sz.nYears, j in 1:sz.nFirms
        if alld_adjust[i, j] == alld[i, j]
            adjust_indicator[i, j] = 1
        end
    end

    return (
        v = allv,
        d = alld,
        a = alla,
        ex = alle,
        y = ally,
        d_adjust = alld_adjust,
        adjust_indicator = adjust_indicator,
        c = allc
    )
end
