function makegrids(ppp::Vector{Float64});
    rho     = ppp[3];
    sigma   = ppp[5];
    rho_e   = ppp[4];
    sigma_e = ppp[6];
    
    #declare the transition matrix and p grid as local and make them
    local trans, trans_p, trans_e, trans_w;

    # Durable Prices
    mew = 0.0; 
    nump = sz.np
    numstd =  sz.nstd
    pg, trans_p = tauchen(mew,sigma,rho,nump,numstd);

    pg = exp.(pg)*1000;

    # Exchange Rate 
    nume = sz.ne
    numstd_e =  sz.nstd_e
    eg, trans_e = tauchen(mew,sigma_e,rho_e,nume,numstd_e);
    eg = exp.(eg);

    #Wage 
    trans_w = [0.8 0.2; 0.2 0.8];
    wg = zeros(2);
    wg[1] = 1000;
    wg[2] = 5000;

    #Kronecker 
    trans = kron(trans_p', kron(trans_e', trans_w'))
    
    #make the durable grids
    dmin = 0.0
    dmax = 50000
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
    amax = 1000000 
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

    outtuple = (p=pg::Vector{Float64},t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64}, e = eg::Vector{Float64}, w = wg::Vector{Float64});
    return outtuple::NamedTuple{(:p,:t,:a,:ap,:d,:dp, :e, :w)};
end

