function simmodel(answ::NamedTuple)
    # Extract necessary objects and grids
    v = answ.v
    pol = answ.pol
    grids = answ.g
    tmat = grids.t
    d_adjust = answ.adjust_result.pol.d
    
    # Initialize the outputs
    allv = zeros(sz.nYears, sz.nFirms)
    alla = zeros(sz.nYears, sz.nFirms)
    alle = zeros(sz.nYears, sz.nFirms)
    ally = zeros(sz.nYears, sz.nFirms)
    alld = zeros(sz.nYears, sz.nFirms)
    alld_adjust = zeros(sz.nYears, sz.nFirms)
    allc = zeros(sz.nYears, sz.nFirms)

    # Precompute cdf for transition matrices
    phatcdf = cumsum(tmat, dims=2)
    
    # Generate aggregate shocks `e` for each year
    # Assume a draw for each year using a random uniform number
    e_unique = rand(1:sz.ne, sz.nYears+1)

    # Preallocate arrays for picke and picky
    all_picke = zeros(Int,sz.nYears+1, sz.nFirms)
    all_picky = zeros(Int, sz.nYears+1, sz.nFirms)

      # Fill in all_picke with the aggregate shock indices
      for ifi in 1:sz.nFirms
        for iti in 1:sz.nYears+1
            # Assign the same aggregate index to all firms for each year
            all_picke[iti, ifi] = e_unique[iti]
        end
    end
    
    # For each firm, determine the initial idiosyncratic shock index randomly
    for ifi in 1:sz.nFirms
        # Initial condition for each firm
        y_draw = globals.draws[1, ifi] .- cumsum(grids.y)
        all_picky[1, ifi] = findfirst(y_draw .< 0.0)
    end
    
    for ifi in 1:sz.nFirms
        for iti in 1:sz.nYears+1
            # Assign the same aggregate shock for each firm
            all_picke[iti, ifi] = e_unique[iti]
            
            # Update for next time period
            if iti <= sz.nYears
                y_draw = globals.draws[iti+1, ifi] .- phatcdf[Int(all_picky[iti,ifi]), :]
                y_draw = y_draw .< 0.0
                combined_index = findfirst(y_draw)

                # Convert the combined index to a specific `y` state index
                if !isnothing(combined_index)
                    all_picky[iti+1, ifi] = mod(combined_index - 1, sz.ny) + 1
                else
                    # Handle case where no valid state found (optional, e.g., set to 1, or handle with better logic)
                    all_picky[iti+1, ifi] = 1
                end


            end
        end
    end
    
    # This is the actual simulation 
    global astart = globals.draws[1, :]
    global dstart = globals.draws[2, :]


    # Perform simulation
    for ifi in 1:sz.nFirms
        # Initialize start positions based on the initial draw
        picka = min(Int(floor(sz.na * astart[ifi])) + 1, sz.na)
        pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
    
        # Using picke and picky from precomputed values
        picke = all_picke[1, ifi]
        picky = all_picky[1, ifi]
    
        vold = v[picke, picky, picka, pickd]
        aold = pol.a[picke, picky, picka, pickd]
        dold = pol.d[picke, picky, picka, pickd]
        d_adjustold = d_adjust[picke, picky, picka, pickd]
        cold = pol.c[picke, picky, picka, pickd]
        
        for iti in 1:sz.nYears
            # Get aggregate and idiosyncratic shocks
            eold = grids.ex[all_picke[iti, ifi]]
            yold = grids.y[all_picky[iti, ifi]]
            
            # Simulate using interpolated policy
            aprime = interpol(eold, yold, aold, dold, grids, pol.a)
            dprime = interpol(eold, yold, aold, dold, grids, pol.d)
            d_adjustprime = interpol(eold, yold, aold, dold, grids, d_adjust)
            cprime = interpol(eold, yold, aold, dold, grids, pol.c)
            vprime = interpol(eold, yold, aold, dold, grids, v)
            
            # Get next period shock indices
            picke = all_picke[iti + 1, ifi]
            picky = all_picky[iti + 1, ifi]
            
            eprime = grids.ex[picke]
            yprime = grids.y[picky]
            
            # Store results
            allv[iti, ifi] = vprime[1]
            alla[iti, ifi] = aprime[1]
            alle[iti, ifi] = eprime
            ally[iti, ifi] = yprime
            alld[iti, ifi] = dprime[1]
            alld_adjust[iti, ifi] = d_adjustprime[1]
            allc[iti, ifi] = cprime[1]
            
            vold = vprime[1]
            aold = aprime[1]
            dold = dprime[1]
            d_adjustold = d_adjustprime[1]
            cold = cprime[1]
        end
    end
    
    # Adjustment indicator
    adjust_indicator = zeros(Int, sz.nYears, sz.nFirms)
    for i in 1:sz.nYears, j in 1:sz.nFirms
        adjust_indicator[i, j] = alld_adjust[i, j] == alld[i, j] ? 1 : 0
    end
    
    return (v=allv, d=alld, a=alla, ex=alle, d_adjust=alld_adjust,
            adjust_indicator=adjust_indicator, c=allc, y=ally, ls=all_picke, picke=all_picke, picky=all_picky)
end