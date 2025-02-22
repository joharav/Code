function makegrids(ppp::Vector{Float64});

    delta   = ppp[2];
    rho_e   = ppp[3];
    sigma_e = ppp[4];
    chi     = ppp[9];
    
    #declare the transition matrix and p grid as local and make them
    local trans;

    # Exchange Rate 
    nume = sz.ne
    numstd_e =  sz.nstd_e
    mew=10.0;
    eg, trans= tauchen(mew,sigma_e,rho_e,nume,numstd_e);
    eg = exp.(eg);
   

    # Durable Grid (State)
    dmin = 0.1
    dmax = 100.0
    if sz.nd == 1
        dg = [0.0]
    else
        dg = collect(range(dmin, stop=dmax, length=sz.nd))  # Durable states
    end

    # Policy Grid (Decision Durable Grid)
    if sz.npd < 2 * sz.nd
        error("npd must be at least 2 * nd to include both adjusted and non-adjusted values.")
    end

    dpg_nonadjust = (1 .- delta .* (1 .- chi)) .* dg  # Apply non-adjustment
    dpg_adjust = collect(range(dmin, stop=dmax, length=sz.nd))  # New durable choices

    dpg = vcat(dpg_nonadjust, dpg_adjust)  # Combine both
    dpg = sort(unique(dpg))  # Remove duplicates and ensure order

    # Define additional grid points to reach npd = 51
    extra_points_needed = sz.npd - length(dpg)

    # Option 1: Add more interpolated points
    if extra_points_needed > 0
        extra_points = collect(range(minimum(dpg)+0.01, stop=maximum(dpg)-0.01, length=extra_points_needed))
        dpg = sort(unique(vcat(dpg, extra_points)))  # Merge and sort
    end
    

    #make the asset grids
    amin = 0.0
    amax = 300.0
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

   
    outtuple = (t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64}, ex = eg::Vector{Float64});
   
    return outtuple::NamedTuple{(:t,:a,:ap,:d,:dp,:ex)};

end
