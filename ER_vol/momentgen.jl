function momentgen(p::Vector{Float64});

    # ============ Run stuff ===================================
    commence = time();
    answ = valfun(p);
    arret = time();
    println("elapse of time in seconds = ",arret-commence)

    if answ.e == 0;
        simdata = simmodel(answ);
        if settings.compstat
            moms, x_values, f_x, h_x, gap_vec = makemoments(simdata,p);
        else
            moms = makemoments(simdata,p);
        end
    else;
        moms = -100.0*ones(sz.nmom,1);
    end;

    return moms::Vector{Float64}; 

    if settings.compstat
        return moms::Vector{Float64}, x_values, f_x, h_x, gap_vec; 
    else
        return moms::Vector{Float64}; 
    end

end