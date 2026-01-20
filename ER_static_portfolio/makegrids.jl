function makegrids(ppp::Vector{Float64})
    # ==========================================================================
    # 4D MODEL GRIDS
    # State: (e, y, w, d) where w = total liquid wealth in pesos
    # Policy: (w', d') with within-period choice of s (dollar share)
    # ==========================================================================
    
    delta = ppp[2]
    rho_e = ppp[3]
    sigma_e = ppp[4]
    rho_y = ppp[14]
    sigma_y = ppp[15]
    pd = ppp[10]
    w_wage = ppp[8]
    tau = ppp[12]
    h = ppp[13]
    chi = ppp[16]

    # -------------------------------------------------------------------------
    # Exchange Rate Grid (discrete Markov)
    # -------------------------------------------------------------------------
    if sz.ne == 2 
        trans_e = [0.90 0.1; 0.85 0.15];
        eg = [1.0, 2.0];
    else
        nume = sz.ne
        numstd_e = sz.nstd_e
        mew = 0.0
        eg, trans_e = tauchen(mew, sigma_e, rho_e, nume, numstd_e)
        eg = exp.(eg)  # ensure E > 0
    end

    # -------------------------------------------------------------------------
    # Income Grid (discrete Markov)
    # -------------------------------------------------------------------------
    if sz.ny == 2 
        trans_y = [0.9 0.1; 0.8 0.2];
        yg = [0.90, 1.1];
    else
        numy = sz.ny
        numstd_y = sz.nstd_y
        mew = 0.0
        yg, trans_y = tauchen(mew, sigma_y, rho_y, numy, numstd_y)
        yg = exp.(yg)
    end
    
    # Joint transition matrix
    trans = kron(trans_e, trans_y)

    # -------------------------------------------------------------------------
    # Total Wealth Grid (replaces separate aa and a grids)
    # w = aa + e*a (total liquid wealth in pesos)
    # -------------------------------------------------------------------------
    income_max = w_wage * h * (1 - tau) * maximum(yg)
    
    # Total wealth grid - should be large enough to capture both asset types
    w_max = max(20.0 * income_max, 10000.0)  
    wg = collect(range(0.0, stop=w_max, length=sz.nw))
    
    # Total wealth policy grid
    wpg = collect(range(0.0, stop=w_max, length=sz.npw))
    
    # -------------------------------------------------------------------------
    # Dollar Share Grid (within-period portfolio choice)
    # s ∈ [0, 1] where s = e*a / w
    # -------------------------------------------------------------------------
    sg = collect(range(0.0, stop=1.0, length=sz.ns))

    # -------------------------------------------------------------------------
    # Durable Grid
    # -------------------------------------------------------------------------
    dmax = max(400.0, 6 * income_max / pd)
    dg = collect(range(0.0, stop=dmax, length=sz.nd))
    
    # Durable policy grid (includes non-adjust continuation points)
    dpg_nonadjust = (1.0 .- delta .* (1.0 .- chi)) .* dg
    dpg_adjust = collect(range(0.0, stop=dmax, length=sz.nd))
    
    dpg = sort(unique(vcat(dpg_nonadjust, dpg_adjust)))
    extra_needed = sz.npd - length(dpg)
    if extra_needed > 0
        extra = collect(range(minimum(dpg)+0.01, stop=maximum(dpg)-0.01, length=extra_needed))
        dpg = sort(unique(vcat(dpg, extra)))
    end
    dpg = dpg[1:min(length(dpg), sz.npd)]  # trim if too many

    if settings.verbose 
        println("=== 4D Model Grids ===")
        println("Income max: ", income_max)
        println("Wealth max: ", w_max)
        println("Durable max: ", dmax)
        println("Total state points: ", sz.ne * sz.ny * sz.nw * sz.nd)
        println("Total policy points per state: ", sz.npw * sz.npd * sz.ns)
    end 

    outtuple = (
        t   = trans::Array{Float64},      # (you can keep it, but stop using it in VFI here)
        ty  = trans_y::Array{Float64},    # <-- ADD THIS
        w   = wg::Vector{Float64},
        wp  = wpg::Vector{Float64},
        s   = sg::Vector{Float64},
        d   = dg::Vector{Float64},
        dp  = dpg::Vector{Float64},
        ex  = eg::Vector{Float64},
        y   = yg::Vector{Float64},
        te  = trans_e::Array{Float64},
    )
    
    return outtuple
end
