function simmodel(answ::NamedTuple);
    #make and extract all of the variables and arrays you need


    v     = answ.v;
    pol   = answ.pol
    grids = answ.g;
    pg    = grids.p
    tmat  = grids.t;


    #initialize the outputs
    allv = zeros(sz.nYears,sz.nFirms);
    alla = zeros(sz.nYears,sz.nFirms);
    allp = zeros(sz.nYears,sz.nFirms);
    alld = zeros(sz.nYears,sz.nFirms);

    #set up the transition matrix
    phatcdf = cumsum(tmat,dims=2);
    
    cdf_wgt = tmat'^100.
    cdf_wgt = cumsum(cdf_wgt[:,Int(floor(sz.np*0.5))+1])

    # ls is the indices of the p shocks in the simulation. stands for "locations" 
    # This sets the starting locations from the unconditional distribution
    ls = zeros(Int, sz.nYears+1, sz.nFirms)
    for ifi in 1:sz.nFirms;
        gap = globals.draws[1,ifi] .- cdf_wgt;
        gap = gap .< 0.0;
        ls[1,ifi] = findfirst(gap);
    end;

    #This is the actual simulation 
    global astart = rand(sz.nFirms);
    global dstart = rand(sz.nFirms);
    for ifi in 1:sz.nFirms;

        #Pick the starting point for the time series 
        picka = min(Int(floor(sz.na*astart[ifi]))+1,sz.na);
        pickd = min(Int(floor(sz.nd*dstart[ifi]))+1,sz.nd);
        pickp = Int(ls[1,ifi]);
        vold =    v[pickp,picka,pickd];
        aold = pol.a[pickp,picka,pickd];
        dold = pol.d[pickp,picka,pickd];
        
        for iti in 1:sz.nYears
            pold = pg[Int(ls[iti,ifi]),1];
            #This updates the simulated variables using simple interpolation
            vprime = interpol(pold,aold,dold,grids,v);  
            aprime = interpol(pold,aold,dold,grids,pol.a); 
            dprime = interpol(pold,aold,dold,grids,pol.d);

            #This updates the shock index using the transition matrix 
            gap = globals.draws[iti+1, ifi] .- phatcdf[Int(ls[iti, ifi]),:];
            gap = gap .< 0.0;
            ls[iti+1,ifi] = findfirst(gap);

            #Update and store
            pprime = pg[Int(ls[iti+1,ifi]),1];
            allv[iti,ifi] =vprime[1]
            alla[iti,ifi] =aprime[1]
            allp[iti,ifi] =pprime
            alld[iti,ifi] =dprime[1]
            vold=vprime[1]
            aold=aprime[1]
            dold=dprime[1]
        end
    end
    outtuple = (v=allv::Array{Float64}, d=alld::Array{Float64}, a=alla::Array{Float64}, p=allp::Array{Float64})
    return outtuple::NamedTuple{(:v,:d,:a,:p)};
end
