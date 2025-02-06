function makegrids(ppp::Vector{Float64});
    rho   = ppp[3];
    sigma = ppp[4];

     #declare the transition matrix and p grid as local and make them
    local trans, p;

    mew = 0.0; 
    nump = sz.np
    numstd =  sz.nstd
    pg, trans = tauchen(mew,sigma,rho,nump,numstd);

    pg = exp.(pg);





    #make the durable grids
    dmin = 0.0
    dmax = 50
    if sz.nd == 1;
        dg = [0.0];
    else;
        dg = collect(range(dmin,stop=dmax,length=sz.nd));
    end;
    if sz.npd == 1;
        dpg = [0.0];
    else;
        dpg = collect(range(dmin,stop=dmax,length=sz.npd));
    end;    


    #make the asset grids
    amin = 0.0
    amax = 100 
    if sz.na == 1;
        ag = [0.0];
    else;
        ag = collect(range(amin,stop=amax,length=sz.na));
    end;
    if sz.npa == 1;
        apg = [0.0];
    else;
        apg = collect(range(amin,stop=amax,length=sz.npa));
    end;    

    outtuple = (p=pg::Vector{Float64},t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64});
    return outtuple::NamedTuple{(:p,:t,:a,:ap,:d,:dp)};
end

