function simmodel(answ::NamedTuple)
    # Make and extract all of the variables and arrays you need
    v           = answ.v
    pol         = answ.pol
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d
    #pol_d         = answ.pol.d

    # Initialize the outputs
    allv        = zeros(sz.nYears, sz.nFirms)
    alla        = zeros(sz.nYears, sz.nFirms)
    alle        = zeros(sz.nYears, sz.nFirms)
    alld        = zeros(sz.nYears, sz.nFirms)
    alld_adjust = zeros(sz.nYears, sz.nFirms)
    allc        = zeros(sz.nYears, sz.nFirms)
    # Set up the transition matrix
    phatcdf = cumsum(tmat, dims=2)
    
    cdf_wgt = tmat'^100
    cdf_wgt = cumsum(cdf_wgt[:,Int(floor(sz.ne*0.5))+1])

    # ls is the indices of the (p, e) shocks in the simulation. stands for "locations" 
    # This sets the starting locations from the unconditional distribution
    ls = zeros(sz.nYears+1,sz.nFirms);
    for ifi in 1:sz.nFirms;
        gap = globals.draws[1,ifi] .- cdf_wgt;
        gap = gap .< 0.0;
        ls[1,ifi] = findfirst(gap);
    end;

    # This is the actual simulation 
    global astart = rand(sz.nFirms)
    global dstart = rand(sz.nFirms)
    
    Threads.@threads for ifi in 1:sz.nFirms
        # Pick the starting point for the time series 
        picka = min(Int(floor(sz.na * astart[ifi])) + 1, sz.na)
        pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
        picke = Int(ls[1, ifi])

        vold =    v[picke,picka,pickd];
        aold = pol.a[picke,picka,pickd];
        dold = pol.d[picke,picka,pickd];
        cold = pol.c[picke,picka,pickd];
        d_adjustold = d_adjust[picke,picka,pickd];
        
        for iti in 1:sz.nYears
            eold = grids.ex[Int(ls[iti,1]),1];

            #This updates the simulated variables using simple interpolation
            vprime = interpol(eold,aold,dold,grids,v);  
            aprime = interpol(eold,aold,dold,grids,pol.a); 
            dprime = interpol(eold,aold,dold,grids,pol.d);
            d_adjustprime = interpol(eold,aold,dold,grids,d_adjust);
            cprime = interpol(eold,aold,dold,grids,pol.c);

            #This updates the shock index using the transition matrix 
            gap = globals.draws[iti+1, ifi] .- phatcdf[Int(ls[iti, ifi]),:];
            gap = gap .< 0.0;
            ls[iti+1,ifi] = findfirst(gap);

            #Update and store
            eprime = grids.ex[Int(ls[iti+1,1]),1];

            allv[iti,ifi]           = vprime[1]
            alla[iti,ifi]           = aprime[1]
            alle[iti,ifi]           = eprime
            alld[iti,ifi]           = dprime[1]
            alld_adjust[iti,ifi]    = d_adjustprime[1]
            allc[iti,ifi]           = cprime[1]


            vold=vprime[1]
            aold=aprime[1]
            dold=dprime[1]
            d_adjustold = d_adjustprime[1]
            cold = cprime[1]
        end
    end

    #check adjustment ratios
    adjust_indicator=zeros(size(alld))
    for i in 1:sz.nYears, j in 1:sz.nFirms
        if alld_adjust[i, j] == alld[i, j]
            adjust_indicator[i, j] = 1
        end
    end

    outtuple = (v=allv::Array{Float64}, d=alld::Array{Float64}, a=alla::Array{Float64}, ex=alle::Array{Float64},d_adjust=alld_adjust::Array{Float64},adjust_indicator=adjust_indicator::Array{Float64},c=allc::Array{Float64})
    return outtuple::NamedTuple{(:v, :d, :a, :ex, :d_adjust, :adjust_indicator, :c)}
end
