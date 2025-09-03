function buildparam(p::Vector{Float64})
    pea = zeros(sz.nop);
    pp = collect(p)
    pea[1:4] = collect(readdlm("../Problemset4/Output/estfil.txt"))  
    pea[5] = pp[1];
 
    return pea::Vector{Float64};
end

# This is a function that is similar to buildparam, but the first four elements of the vector and keeps the fifth at its original value. 
function buildparam(p::Float64,sub_pea::Vector{Float64})
    pp = collect(p)
    pea = zeros(sz.nop);
    pea[1:4] = sub_pea[1:4] 
    pea[5] = pp[1];
 
    return pea::Vector{Float64};
end


# ============== The fcn function has two methods. One returns moments and basically wraps a wrapper. 
# ===============The other returns a GMM objective function. 

function fcn(p::Vector{Float64})
    pea = buildparam(p);
    moms = momentgen(pea);
    moms = moms[sz.pick]
    return moms::Vector{Float64}
end 

function subfcn(subp::Vector{Float64})
    pp = collect(readdlm("Output/estfil.txt"))
    pp = pp[1]  
    pea = buildparam(pp,subp);
    moms = momentgen(pea);
    moms = moms[sz.pick]
    return moms::Vector{Float64}
end 


function fcn(p::Vector{Float64},fopt::Float64);

    pea = buildparam(p)
    
    simmoms     = collect(momentgen(pea));

    if simmoms[1] != -100.0
        datamoms    = collect(readdlm("../Problemset3/moments.txt")); 
        sw          = collect(readdlm("../Problemset3/cluster_covariance_matrix.txt"));
        momname     = collect(readdlm("../Problemset3/momnames.txt"));
        pname       = collect(readdlm("../Problemset3/pnames.txt"));

        ch = sz.pick;#.==ones(sz.nmom)); println(ch); println(" "); println(sw;)
        sw = sw[ch,ch]
        
        if settings.complicated
            w = I(size(sw,1))
        else
            w = inv(sw);
        end
        datamoms = datamoms[ch]
        simmoms  = simmoms[ch]; 
        momdiff  = datamoms .- simmoms;
        momname  = momname[ch]

        bigQ     = momdiff'*w*momdiff;
        bigQ    = bigQ[1];
        #println(bigQ,typeof(bigQ),size(bigQ))
        #println(simmoms,typeof(simmoms),size(simmoms))
        #println(datamoms,typeof(datamoms),size(datamoms))
        #println(sw[:,1],typeof(sw[:,1]),size(sw[:,1]))
 
        if bigQ < fopt;
            filname = "Output/progress.txt"
            io = open(filname,"w")
            @printf(io," data and simulated moments\n")
            for jj in 1:size(momdiff,1)
                @printf(io,"%s",momname[jj]);
                @printf(io,"%16.8f",datamoms[jj])  
                @printf(io,"%16.8f\n",simmoms[jj])  
            end
            @printf(io," \n")
            @printf(io," parameter value\n")
            for jj in eachindex(p);
                @printf(io,"%s","lambda"); 
                @printf(io,"%16.8f\n",p[jj])
            end
            @printf(io," \n")
            @printf(io,"GMM function value =   %16.8f\n",bigQ)
            close(io)
        end;
    else; 
        bigQ = sz.maxfunc
    end

    return bigQ::Float64
end

function smmstats(p::Vector{Float64})
    pea = buildparam(p)

    simmoms     = collect(momentgen(pea));
    datamoms    = collect(readdlm("../Problemset3/moments.txt")); 
    sw          = collect(readdlm("../Problemset3/cluster_covariance_matrix.txt"));
    gvkey       = collect(readdlm("../Problemset3/gvkeys.txt"));
    infl_fcn    = collect(readdlm("../Problemset3/influence_functions.txt")); #data calculations. moments 
    infl_param  = collect(readdlm("../Problemset4/Output/influence_functions.txt")); #parameter influence functions

    infl_fcn    = infl_fcn[sz.pick,:]

    #println(size(sw))
    ch = sz.pick;
    sw = sw[ch,ch]
    datamoms = datamoms[ch]
    simmoms  = simmoms[ch]

    infl_fcn_for_params = infl_fcn .+ datamoms .- simmoms

    samplesize  = size(gvkey,1);

    if settings.complicated
        w = I(size(sw,1))
    else
        w = inv(sw);
    end

    gee = grad(p, sz.nmom, sz.noestp); println(gee)
  
    

    gwg   = gee'*w*gee;                                              
    igwg = inv(gwg);                                                   

    # This part does the plug-in correction. 

    sub_pea = vec(collect(readdlm("../Problemset4/Output/estfil.txt")))  
    sub_pea = sub_pea[1:4]



    jee = grad(sub_pea, sz.nmom, size(sub_pea,1));  println(jee)

    jee = jee'
    println(typeof(jee),size(jee))
    println(typeof(infl_param),size(infl_param))
    println(typeof(infl_fcn),size(infl_fcn))

    ugly_infl = infl_fcn' - infl_param*jee

    sw = ugly_infl'*ugly_infl/samplesize/samplesize;

    wsw   = w*sw*w;                                                   
    gwswg = gee'*wsw*gee;                                               

    println("sw size is ",size(sw))

     #This part makes a new influence function by 

    if settings.complicated                                              
       vc    =  igwg*gwswg*igwg;                                      
    else                                                              
       vc    =  igwg;                                                 
    end                                                             
    vc    =  igwg*gwswg*igwg;  #This is kludgy. It's just for a problem set though. 
    
    nsim = float((sz.nYears-sz.burnin)*sz.nFirms)/float(samplesize)
    vc    = (1.0 + 1.0/nsim)*vc # for code where you don't divide by n^2 in making the weight matrix, you divide by the sample size here

    if sz.noestp == 1;
        standarderror = sqrt(vc)
    else
        standarderror = sqrt.(diag(vc))
    end

 println("size of gee is ",size(gee))
 println("size of igwg is ",size(igwg))


    gigwgg = gee*igwg*gee';  
    eye = float(collect(I(size(w,1))))
    eyegg = eye - gigwgg*w;
    println("eyegg size is, ", size(eyegg)) 
    vpe = eyegg*sw*eyegg'; 
    vpe = (1.0 + 1.0/nsim)*vpe # for code where you don't divide by n^2 in making the weight matrix, you divide by the sample size here
    ivpe = pinv(vpe);  # Moore Penrose
    
    momdiff = datamoms-simmoms
   
    if settings.complicated 
        jtest = momdiff'*ivpe*momdiff;  
    else
        jtest = momdiff'*w*momdiff; # for code where you don't divide by n^2 in making the weight matrix, you muliply by the sample size here
    end
    jtest = momdiff'*ivpe*momdiff;  #This is also kludgy. It's just for a problem set though.

    genshap =  igwg*gee'*w; 

    for jj = 1:sz.nmom
        for ii = 1:sz.noestp
            genshap[ii,jj] = genshap[ii,jj] * sqrt(sw[ii,ii] / samplesize)
        end
    end
    infl_parameters = genshap*infl_fcn_for_params
    
   outtuple = (datamoms=datamoms, simmoms=simmoms, standarderror=standarderror, vpe=vpe, gee=gee, jtest=jtest, genshap=genshap, infl=infl_parameters)
   return outtuple 

end

function print_smm_results(smm_results,p)

    datamoms        = smm_results.datamoms
    simmoms         = smm_results.simmoms
    standarderror   = smm_results.standarderror
    vpe             = smm_results.vpe
    gee             = smm_results.gee
    genshap         = smm_results.genshap
    jtest           = smm_results.jtest[1]
    infl            = smm_results.infl

    momname     = collect(readdlm("../Problemset3/momnames.txt"));
    pname       = collect(readdlm("../Problemset3//pnames.txt"));
    momname    = momname[sz.pick]
    
    if settings.complicated;
        filname = "Output/identity_results"
    else
        filname = "Output/results.txt"
    end

    io = open(filname,"w")
    @printf(io,"Final SMM results in LaTeX format \n")
    @printf(io," \n")

    @printf(io,"Parameter estimates and standard errors \n")
    @printf(io," \n")    
    
    for jj in eachindex(p); @printf(io,"&  %s   ",pname[jj]); end;   @printf(io," \n") 
    for jj in eachindex(p); @printf(io,"&  %12.4f  ",p[jj]);  end;   @printf(io," \n")  
    for jj in eachindex(p); @printf(io,"&( %12.4f )",standarderror[jj]); end;   @printf(io," \n")  
    @printf(io," \n")  

    @printf(io,"Data and Simulated Moments and t-stats \n")
    @printf(io," \n")    
    for jj in eachindex(datamoms);
        tstat = (datamoms[jj]-simmoms[jj])/sqrt(vpe[jj,jj]);
        @printf(io,"%s & %12.4f & %12.4f & %12.4f \\\\ \n",momname[jj],datamoms[jj],simmoms[jj],tstat); 
    end;   
    @printf(io," \n")  
        
    @printf(io,"Overidentification test \n")
    @printf(io," \n")    

    @printf(io,"GMM J-test %12.4f",jtest);
    @printf(io," \n")  
    
    @printf(io," \n")  
    
    @printf(io,"Andrews Gentzgow Shapiro \n")
    @printf(io," \n")    


    @printf(io,"           ")
    for jj in 1:sz.noestp; @printf(io,"&   %s   ",pname[jj]); end;   @printf(io," \n") 
    for ii in 1:sz.nmom
    @printf(io,"%s",momname[ii])
    for jj in eachindex(p);   ;@printf(io,"&  %12.4f  ",genshap[jj,ii]);  end;   @printf(io," \n")  
    end;
 
    @printf(io," \n")
    @printf(io,"Jacobian matrix \n")
    @printf(io," \n")    
    for ii in 1:sz.nmom; 
    for jj in 1:sz.noestp; @printf(io,"&  %12.4f  ",gee[ii,jj]);  end;   @printf(io," \n")  
    end;
    close(io)

    filname = "Output/influence_functions.txt"
    io = open(filname,"w")
    for ii in 1:size(infl,2)
    for jj in 1:sz.noestp; @printf(io,"  %36.16f  ",infl[jj,ii]);  end;   @printf(io," \n")
    end;
    close(io)
end;


function grad(x0, n, k)
    g = zeros(n, k)
    
    ax0 = abs.(x0)
    dax0 = sign.(x0)
    dh = (1e-1) .* ax0 .* dax0
    xdup = x0 + dh
    xddw = x0 - dh
    
    argup = repeat(x0', k, 1)
    argdw = repeat(x0', k, 1)
    
    for i in 1:k
        argup[i, i] = xdup[i]
        argdw[i, i] = xddw[i]
    end
    
    g_up = zeros(n, k)
    g_dw = zeros(n, k)
    
    for i in 1:k
        g_up[:, i] = fcn(argup[:, i])
        g_dw[:, i] = fcn(argdw[:, i])
    end
    
    for i in 1:k
        g[:, i] = (g_up[:, i] - g_dw[:, i]) / (2.0 * dh[i])
    end
    
    return g
end