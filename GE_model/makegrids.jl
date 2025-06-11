function makegrids(ppp::Vector{Float64})
    delta = ppp[2]
    rho_e = ppp[3]
    sigma_e = ppp[4]
    chi = ppp[5]


    # Exchange Rate (as depreciation shocks)
    if sz.ne == 1
        eg = [1.0]
        trans_e = [1.0]
    else
        nume = sz.ne
        numstd_e = sz.nstd_e
        mew = 0.0  # mean log depreciation = 0 → implies median e = 1
        eg_log, trans_e = tauchen(mew, sigma_e, rho_e, nume, numstd_e)

        # eg_log now contains log-depreciation rates: log(e_t)
        # Convert to levels: gross depreciation rate
        eg = exp.(eg_log)

        # Ensure middle point is exactly 1 (normalize)
        mid_idx = Int(ceil(nume / 2))
        eg ./= eg[mid_idx]
    end
    
    
    if sz.nz == 3 
        trans_z = [0.5 0.4 0.1; 0.3 0.5 0.2 ; 0.1 0.5 0.4];
        zg = zeros(3);
        zg[1] = 0.90;
        zg[2] = 1.0;
        zg[3] = 1.10;
    else

        numz = sz.nz
        numstd_z = sz.nstd_z
        mew = 0.0
        zg, trans_z = tauchen(mew, sigma_z, rho_z, numz, numstd_z)
        zg=exp.(zg)

    end
        #Kronecker 
        trans = kron(trans_z, trans_e)



    ### Non-Adjust Case: Asset and Durable Grids

    # Asset Grid (Current Asset Holdings)
    a_min = 0.0
    a_max = 50
    ag = collect(range(a_min, stop=a_max, length=sz.na))
    apg = collect(range(a_min, stop=a_max, length=sz.npa))

    # Durable Grid (State)
    dmin = 0.0
    dmax = 60.0
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



    outtuple = (t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64}, ex = eg::Vector{Float64}, z = zg::Vector{Float64});
   
    return outtuple::NamedTuple{(:t,:a,:ap,:d,:dp,:ex,:z)};
end

