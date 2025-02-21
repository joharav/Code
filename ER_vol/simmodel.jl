function simmodel(answ::NamedTuple)
    # Make and extract all of the variables and arrays you need
    v     = answ.v
    pol   = answ.pol
    grids = answ.g
    tmat  = grids.t

    # Initialize the outputs
    allv = zeros(sz.nYears, sz.nFirms)
    alla = zeros(sz.nYears, sz.nFirms)
    allp = zeros(sz.nYears, sz.nFirms)
    alle = zeros(sz.nYears, sz.nFirms)
    alld = zeros(sz.nYears, sz.nFirms)

    # Set up the transition matrix
    phatcdf = cumsum(tmat, dims=2)
    
    cdf_wgt = tmat'^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.np * sz.ne * 0.5)) + 1])

    # ls is the indices of the (p, e) shocks in the simulation. stands for "locations" 
    # This sets the starting locations from the unconditional distribution
    ls = zeros(Int, sz.nYears + 1, sz.nFirms)
    Threads.@threads for ifi in 1:sz.nFirms
        gap = globals.draws[1, ifi] .- cdf_wgt
        gap = gap .< 0.0
        ls[1, ifi] = findfirst(gap)
    end

    # This is the actual simulation 
    global astart = rand(sz.nFirms)
    global dstart = rand(sz.nFirms)
    Threads.@threads for ifi in 1:sz.nFirms
        # Pick the starting point for the time series 
        picka = min(Int(floor(sz.na * astart[ifi])) + 1, sz.na)
        pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
        pick_combined = Int(ls[1, ifi])
        picke = Int(floor((pick_combined - 1) / (sz.np)) % sz.ne) + 1
        pickp = Int((pick_combined - 1) % sz.np) + 1

        # Store the initial values
        alla[1, ifi] = pol.a[pickp, picke, picka, pickd]
        alld[1, ifi] = pol.d[pickp, picke, picka, pickd]
        allv[1, ifi] = v[pickp, picke, picka, pickd]
        allp[1, ifi] = pickp
        alle[1, ifi] = picke

        # Simulate over the years
        for t in 2:sz.nYears
            # Draw the next state
            gap = globals.draws[t, ifi] .- phatcdf[pick_combined, :]
            gap = gap .< 0.0
            pick_combined = findfirst(gap)

            # Update the state variables
            picke = Int(floor((pick_combined - 1) / (sz.np)) % sz.ne) + 1
            pickp = Int((pick_combined - 1) % sz.np) + 1

            # This updates the simulated variables using simple interpolation
            vprime = interpol(Float64(pickp), Float64(picke), Float64(picka), Float64(pickd), grids, v)
            aprime = interpol(Float64(pickp), Float64(picke), Float64(picka), Float64(pickd), grids, pol.a)
            dprime = interpol(Float64(pickp), Float64(picke), Float64(picka), Float64(pickd), grids, pol.d)

            # Store the values
            alla[t, ifi] = aprime
            alld[t, ifi] = dprime
            allv[t, ifi] = vprime
            allp[t, ifi] = pickp
            alle[t, ifi] = picke
        end
    end

    outtuple = (v=allv::Array{Float64}, d=alld::Array{Float64}, a=alla::Array{Float64}, p=allp::Array{Float64}, ex=alle::Array{Float64})
    return outtuple::NamedTuple{(:v, :d, :a, :p, :ex)}
end
