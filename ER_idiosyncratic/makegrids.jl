function makegrids(ppp::Vector{Float64})
    delta = ppp[2]
    rho_e = ppp[3]
    sigma_e = ppp[4]
    chi = ppp[5]
    rho_y = ppp[14]
    sigma_y = ppp[15]


    # Exchange Rate
    if sz.ne == 3 
        # Exchange Rate
        trans_e = [0.5 0.4 0.1; 0.3 0.5 0.2 ; 0.1 0.5 0.4];
        eg = zeros(3);
        eg[1] = 0.5;
        eg[2] = 1.0;
        eg[3] = 2.0;
    else
        nume = sz.ne
        numstd_e = sz.nstd_e
        mew = 3.0
        eg, trans_e = tauchen(mew, sigma_e, rho_e, nume, numstd_e)
            # Rescale to (0,2]
        eg = 0.01 .+ (1.99) .* (eg .- minimum(eg)) ./ (maximum(eg) - minimum(eg))
    end



    
    
    # Idiosyncratic income

    if sz.ny == 3 
        trans_y = [0.5 0.4 0.1; 0.3 0.5 0.2 ; 0.1 0.5 0.4];
        yg = zeros(3);
        yg[1] = 0.90;
        yg[2] = 1.0;
        yg[3] = 1.10;
    else

        numy = sz.ny
        numstd_y = sz.nstd_y
        mew = 0.0
        yg, trans_y = tauchen(mew, sigma_y, rho_y, numy, numstd_y)
        yg=exp.(yg)

    end
        #Kronecker 
        trans = kron(trans_e, trans_y)

    ### Non-Adjust Case: Asset and Durable Grids

    # Asset Grid (Current Asset Holdings)
    a_min = 0.0
    a_max = 80
    ag = collect(range(a_min, stop=a_max, length=sz.na))
    apg = collect(range(a_min, stop=a_max, length=sz.npa))

    # Durable Grid (State)
    dmin = 0.0
    dmax = 90.0
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

    if extra_points_needed > 0
    extra_points = collect(range(minimum(dpg)+0.01, stop=maximum(dpg)-0.01, length=extra_points_needed))
    dpg = sort(unique(vcat(dpg, extra_points)))  # Merge and sort
    end



    outtuple = (t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64}, ex = eg::Vector{Float64}, y = yg::Vector{Float64});
   
    return outtuple::NamedTuple{(:t,:a,:ap,:d,:dp,:ex,:y)};
end

