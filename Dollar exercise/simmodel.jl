function simmodel(answ::NamedTuple);
    #make and extract all of the variables and arrays you need

    v     = answ.v;
    pol   = answ.pol
    grids = answ.g;
    pg    = grids.p
    eg    = grids.e
    wg    = grids.w
    tmat  = grids.t;

    #initialize the outputs
    allv = zeros(sz.nYears,sz.nFirms);
    alla = zeros(sz.nYears,sz.nFirms);
    allp = zeros(sz.nYears,sz.nFirms);
    alle = zeros(sz.nYears,sz.nFirms);
    allw = zeros(sz.nYears,sz.nFirms);
    alld = zeros(sz.nYears,sz.nFirms);

    #set up the transition matrix
    phatcdf = cumsum(tmat,dims=2);
    
    cdf_wgt = tmat'^100.
    cdf_wgt = cumsum(cdf_wgt[:,Int(floor(sz.np*sz.ne*sz.nw*0.5))+1])

    # ls is the indices of the p, e, w shocks in the simulation. stands for "locations" 
    # This sets the starting locations from the unconditional distribution
    ls = zeros(Int, sz.nYears+1, sz.nFirms)
    Threads.@threads for ifi in 1:sz.nFirms;
        gap = globals.draws[1,ifi] .- cdf_wgt;
        gap = gap .< 0.0;
        ls[1,ifi] = findfirst(gap);
    end;

    # This is the actual simulation 
    global astart = rand(sz.nFirms);
    global dstart = rand(sz.nFirms);
    Threads.@threads for ifi in 1:sz.nFirms;

        # Pick the starting point for the time series 
        picka = min(Int(floor(sz.na*astart[ifi]))+1,sz.na);
        pickd = min(Int(floor(sz.nd*dstart[ifi]))+1,sz.nd);
        pickp = Int(ls[1,ifi]);
        picke = Int(floor((pickp - 1) / (sz.np * sz.nw)) % sz.ne) + 1
        pickw = Int(floor((pickp - 1) / sz.np) % sz.nw) + 1
        pickp = Int((pickp - 1) % sz.np) + 1
        vold = v[pickp,picke,pickw,picka,pickd];
        aold = pol.a[pickp,picke,pickw,picka,pickd];
        dold = pol.d[pickp,picke,pickw,picka,pickd];
        
        Threads.@threads for iti in 1:sz.nYears
            pold = pg[pickp];
            eold = eg[picke];
            wold = wg[pickw];
            # This updates the simulated variables using simple interpolation
            vprime = interpol(pold, eold, wold, aold, dold, grids, v);  
            aprime = interpol(pold, eold, wold, aold, dold, grids, pol.a); 
            dprime = interpol(pold, eold, wold, aold, dold, grids, pol.d);

            #This updates the shock index using the transition matrix 
            gap = globals.draws[iti+1, ifi] .- phatcdf[Int(ls[iti, ifi]),:];
            gap = gap .< 0.0;
            ls[iti+1,ifi] = findfirst(gap);

            # Update and store
            pprime = pg[Int(ls[iti+1,ifi]), 1]
            eprime = eg[Int(floor((ls[iti+1,ifi] - 1) / (sz.np * sz.nw)) % sz.ne) + 1]
            wprime = wg[Int(floor((ls[iti+1,ifi] - 1) / sz.np) % sz.nw) + 1]
            allv[iti,ifi] = vprime[1]
            alla[iti,ifi] = aprime[1]
            allp[iti,ifi] = pprime
            alle[iti,ifi] = eprime
            allw[iti,ifi] = wprime
            alld[iti,ifi] = dprime[1]
            vold = vprime[1]
            aold = aprime[1]
            dold = dprime[1]
        end
    end
    outtuple = (v=allv::Array{Float64}, d=alld::Array{Float64}, a=alla::Array{Float64}, p=allp::Array{Float64}, e=alle::Array{Float64}, w=allw::Array{Float64})
    return outtuple::NamedTuple{(:v, :d, :a, :p, :e, :w)}
end
