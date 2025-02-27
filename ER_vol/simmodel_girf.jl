#Generalized Impulse Response Function

function simmodel_girf(answ::NamedTuple, T_shock::Int, shock_size::Float64)
    v = answ.v
    pol = answ.pol
    grids = answ.g
    tmat = grids.t
    d_adjust = answ.adjust_result.pol.d

    allv_shock        = zeros(sz.nYears, sz.nFirms)
    alla_shock        = zeros(sz.nYears, sz.nFirms)
    alle_shock        = zeros(sz.nYears, sz.nFirms)
    alld_shock        = zeros(sz.nYears, sz.nFirms)
    alld_adjust_shock = zeros(sz.nYears, sz.nFirms)
    
    phatcdf = cumsum(tmat, dims=2)
    cdf_wgt = tmat'^100
    cdf_wgt = cumsum(cdf_wgt[:,Int(floor(sz.ne*0.5))+1])

    ls = zeros(sz.nYears+1,sz.nFirms)
    for ifi in 1:sz.nFirms
        gap = globals.draws[1,ifi] .- cdf_wgt
        gap = gap .< 0.0
        ls[1,ifi] = findfirst(gap)
    end

        global astart = rand(sz.nFirms)
        global dstart = rand(sz.nFirms)

        Threads.@threads for ifi in 1:sz.nFirms
            picka = min(Int(floor(sz.na * astart[ifi])) + 1, sz.na)
            pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
            picke = Int(ls[1, ifi])

            vold = v[picke,picka,pickd]
            aold = pol.a[picke,picka,pickd];
            dold = pol.d[picke,picka,pickd];
            d_adjustold = d_adjust[picke,picka,pickd];
            
            for iti in 1:sz.nYears
                eold = grids.ex[Int(ls[iti,1]),1]
                
                if iti == T_shock
                    eold *= (1 + shock_size)  # Apply the exchange rate shock
                end

                #This updates the simulated variables using simple interpolation
                vprime = interpol(eold,aold,dold,grids,v);  
                aprime = interpol(eold,aold,dold,grids,pol.a); 
                dprime = interpol(eold,aold,dold,grids,pol.d);
                d_adjustprime = interpol(eold,aold,dold,grids,d_adjust);           

                gap = globals.draws[iti+1, ifi] .- phatcdf[Int(ls[iti, ifi]),:]
                gap = gap .< 0.0
                ls[iti+1,ifi] = findfirst(gap)
                
                #Update and store
                eprime = grids.ex[Int(ls[iti+1,1]),1];

                allv_shock[iti,ifi]           = vprime[1]
                alla_shock[iti,ifi]           = aprime[1]
                alle_shock[iti,ifi]           = eprime
                alld_shock[iti,ifi]           = dprime[1]
                alld_adjust_shock[iti,ifi]    = d_adjustprime[1]


                vold=vprime[1]
                aold=aprime[1]
                dold=dprime[1]
                d_adjustold = d_adjustprime[1]

                end
        end
    

    
    adjust_indicator_shock=zeros(size(alld_adjust_shock))
    for i in 1:sz.nYears, j in 1:sz.nFirms
        if alld_adjust_shock[i, j] == alld_adjust_shock[i, j]
            adjust_indicator_shock[i, j] = 1
        end
    end
    
    outtuple = (v=allv_shock::Array{Float64}, d=alld_shock::Array{Float64}, a=alla_shock::Array{Float64}, ex=alle_shock::Array{Float64},d_adjust_shock=alld_adjust_shock::Array{Float64},adjust_indicator=adjust_indicator_shock::Array{Float64})
    return outtuple::NamedTuple{(:v, :d, :a, :ex, :d_adjust, :adjust_indicator)}
end
