function makegrids(ppp::Vector{Float64})
    delta = ppp[2]
    rho_e = ppp[3]
    sigma_e = ppp[4]
    rho_y = ppp[14]
    sigma_y = ppp[15]
    pd = ppp[10]
    w = ppp[8]
    tau = ppp[12]
    h = ppp[13]
    chi = ppp[16]


    # Exchange Rate
    if sz.ne == 2 
        # Exchange Rate
        trans_e = [0.90 0.1; 0.85 0.15];
        eg = zeros(2);
        eg[1] = 1.0;
        eg[2] = 2.0;
    else
        nume = sz.ne
        numstd_e = sz.nstd_e
        mew = 3.0
        eg, trans_e = tauchen(mew, sigma_e, rho_e, nume, numstd_e)
        eg = exp.(eg)  # ensure E > 0
        
    end




    # Idiosyncratic income

    if sz.ny == 2 
        trans_y = [0.9 0.1; 0.8 0.2];
        yg = zeros(2);
        yg[1] = 0.90;
        yg[2] = 1.1;
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

    income_max = w * h * (1 - tau) * maximum(yg)
    
    # Local Asset Grid
    aa_max = max(10.0 * income_max, 5000.0)  # increased from 800
    aag = collect(range(0.0, stop=aa_max, length=sz.na))
    aapg = collect(range(0.0, stop=aa_max, length=sz.npa))



   # Asset Grid
   a_max = max(10.0 * income_max, 5000.0)  # increased from 800
   ag = collect(range(0.0, stop=a_max, length=sz.na))
   apg = collect(range(0.0, stop=a_max, length=sz.npa))

    # Durable Grid
    dmax = max(400.0, 6 * income_max / pd)
    dg = collect(range(0.0, stop=dmax, length=sz.nd))
    # Durable Policy Grid (combine adjust and non-adjust)
    dpg_nonadjust = (1 .- delta .* (1 .- chi)) .* dg
    dpg_adjust = collect(range(0.0, stop=dmax, length=sz.nd))


    dpg = sort(unique(vcat(dpg_nonadjust, dpg_adjust)))
    extra_needed = sz.npd - length(dpg)
    if extra_needed > 0
        extra = collect(range(minimum(dpg)+0.01, stop=maximum(dpg)-0.01, length=extra_needed))
        dpg = sort(unique(vcat(dpg, extra)))
    end

    if settings.verbose 
        println("Income max: ", income_max)
        println("a_max: ", a_max)
        println("d_max: ", dmax)
    end 


    outtuple = (t=trans::Array{Float64},a=ag::Vector{Float64},ap=apg::Vector{Float64},d=dg::Vector{Float64},dp = dpg::Vector{Float64}, ex = eg::Vector{Float64}, y = yg::Vector{Float64}, aa = aag::Vector{Float64}, aap = aapg::Vector{Float64})
    return outtuple::NamedTuple{(:t,:a,:ap,:d,:dp,:ex,:y,:aa,:aap)}
end


