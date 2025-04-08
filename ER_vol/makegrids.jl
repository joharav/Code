function makegrids(ppp::Vector{Float64})
    delta = ppp[2]
    rho_e = ppp[3]
    sigma_e = ppp[4]
    chi = ppp[5]


    # Exchange Rate
    if sz.ne==1
        eg=[1.0]
        trans = [1.0]
    else
        nume = sz.ne
        numstd_e = sz.nstd_e
        mew = 3.0
        eg, trans = tauchen(mew, sigma_e, rho_e, nume, numstd_e)

        # Rescale to (0,2]
        eg = 0.01 .+ (1.99) .* (eg .- minimum(eg)) ./ (maximum(eg) - minimum(eg))
    end
    
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



    outtuple = (t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64}, ex = eg::Vector{Float64});
   
    return outtuple::NamedTuple{(:t,:a,:ap,:d,:dp,:ex)};
end

