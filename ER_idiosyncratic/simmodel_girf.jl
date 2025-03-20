#Generalized Impulse Response Function

function simmodel_girf(answ::NamedTuple, T_shock::Int)
    v           = answ.v
    pol         = answ.pol
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d

    allv_shock        = zeros(sz.nYears, sz.nFirms)
    alla_shock        = zeros(sz.nYears, sz.nFirms)
    alle_shock        = zeros(sz.nYears, sz.nFirms)
    ally_shock        = zeros(sz.nYears, sz.nFirms)
    alld_shock        = zeros(sz.nYears, sz.nFirms)
    alld_adjust_shock = zeros(sz.nYears, sz.nFirms)
    allc_shock        = zeros(sz.nYears, sz.nFirms)
    
    phatcdf = cumsum(tmat, dims=2)
    cdf_wgt = tmat'^100
    cdf_wgt = cumsum(cdf_wgt[:,Int(floor(sz.ne*sz.ny*0.5))+1])

    ls = zeros(sz.nYears+1,sz.nFirms)
    for ifi in 1:sz.nFirms
        gap = globals.draws[1,ifi] .- cdf_wgt
        gap = gap .< 0.0
        ls[1,ifi] = findfirst(gap)
    end

    for ifi in 1:sz.nFirms
        for iti in 1:sz.nYears

        #This updates the shock index using the transition matrix 
        gap = globals.draws[iti+1, ifi] .- phatcdf[Int(ls[iti, ifi]),:];
        gap = gap .< 0.0;
        ls[iti+1,ifi] = findfirst(gap);
        end
    end 

    # Initialize matrices to store picke and picky values
    all_picke = zeros(Int, sz.nYears+1, sz.nFirms)
    all_picky = zeros(Int, sz.nYears+1, sz.nFirms)

    # Compute picke and picky for each value in ls
    for ifi in 1:sz.nFirms
        for iti in 1:sz.nYears+1
            pick_combined = Int(ls[iti, ifi])  # Extract the combined index

            # Compute the indices
            all_picke[iti, ifi] = div(pick_combined - 1, sz.ny) + 1
            all_picky[iti, ifi] = mod(pick_combined - 1, sz.ny) + 1
        end
    end



    global astart = globals.draws[1, :]
    global dstart = globals.draws[2, :]

        Threads.@threads for ifi in 1:sz.nFirms
            picka = min(Int(floor(sz.na * astart[ifi])) + 1, sz.na)
            pickd = min(Int(floor(sz.nd * dstart[ifi])) + 1, sz.nd)
    
            picke = all_picke[1, ifi]
            picky = all_picky[1, ifi]

            vold            = v[picke,picky,picka,pickd]
            aold            = pol.a[picke,picky,picka,pickd];
            dold            = pol.d[picke,picky,picka,pickd];
            d_adjustold     = d_adjust[picke,picky,picka,pickd];
            cold            = pol.c[picke,picky,picka,pickd];
            
            for iti in 1:sz.nYears
                eold = grids.ex[picke];
                yold = grids.y[picky];

                if iti == T_shock
                    eold = grids.ex[sz.ne]  # Apply the exchange rate shock
                end

                #This updates the simulated variables using simple interpolation
                aprime              = interpol(eold,yold,aold,dold,grids,pol.a); 
                dprime              = interpol(eold,yold,aold,dold,grids,pol.d);
                d_adjustprime       = interpol(eold,yold,aold,dold,grids,d_adjust);     
                cprime              = interpol(eold,yold,aold,dold,grids,pol.c);      
                vprime              = interpol(eold,yold,aold,dold,grids,v);  

                
                #Update and store
                combinedprime = Int(ls[iti+1,ifi])
                
                # Corrected extraction of e' and y' using the same logic as before
                picke = all_picke[iti+1, ifi]
                picky = all_picky[iti+1, ifi]

                eprime = grids.ex[picke]  # **Updated extraction of e'**
                yprime = grids.y[picky]  # **Updated extraction of y'**

                allv_shock[iti,ifi]           = vprime[1]
                alla_shock[iti,ifi]           = aprime[1]
                alle_shock[iti,ifi]           = eprime
                ally_shock[iti,ifi]           = yprime
                alld_shock[iti,ifi]           = dprime[1]
                alld_adjust_shock[iti,ifi]    = d_adjustprime[1]
                allc_shock[iti,ifi]           = cprime[1]

                if iti == T_shock
                    alle_shock[iti,ifi]           = eold # Apply the exchange rate shock
                end

                vold=vprime[1]
                aold=aprime[1]
                dold=dprime[1]
                d_adjustold = d_adjustprime[1]
                cold=cprime[1]

                end
        end
    
    #check adjustment ratios
    adjust_indicator=zeros(size(alld_shock))
    for i in 1:sz.nYears, j in 1:sz.nFirms
        if alld_adjust_shock[i, j] == alld_shock[i, j]
            adjust_indicator[i, j] = 1
        end
    end
    
    outtuple = (v=allv_shock::Array{Float64}, d=alld_shock::Array{Float64}, a=alla_shock::Array{Float64}, ex=alle_shock::Array{Float64},d_adjust=alld_adjust_shock::Array{Float64},adjust_indicator=adjust_indicator::Array{Float64}, c=allc_shock::Array{Float64}, y=ally_shock::Array{Float64})
    return outtuple::NamedTuple{(:v, :d, :a, :ex, :d_adjust, :adjust_indicator, :c, :y)}
end
