function momentgen(p::Vector{Float64});

    # ============ Run stuff ===================================
    commence = time();
    answ = valfun(p);
    arret = time();
    println("elapse of time in seconds = ",arret-commence)
    adjust_result=answ.adjust_result;
    noadjust_result=answ.noadjust_result;

    if answ.e == 0;
        simdata = simmodel(answ);
        moms = makemoments(simdata,p,adjust_result,noadjust_result);
    else;
        moms = -100.0*ones(sz.nmom,1);
    end;

    return moms::Vector{Float64}; 
end