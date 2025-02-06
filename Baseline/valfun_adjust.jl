function valfun_adjust(pea::Vector{Float64}) 
    beta = pea[1]

    #initialize stuff
    v        = zeros(sz.np,sz.na,sz.nd);     #value function
    vnew     = zeros(sz.np,sz.na,sz.nd);     #new value function 
    mew      = zeros(sz.np,sz.na,sz.nd);     #stationary distribution

    #initialize and instantiate the integer policy function
    gidx = dtp.Ipol(Int.(zeros(sz.np, sz.na, sz.nd)), Int.(zeros(sz.np, sz.na, sz.nd)))
    pol = dtp.Pol(zeros(sz.np, sz.na, sz.nd), zeros(sz.np, sz.na, sz.nd))

    #initialize and instantiate the old integer policy function
    old = dtp.Ipol(Int.(zeros(sz.np, sz.na, sz.nd)), Int.(zeros(sz.np, sz.na, sz.nd)))  
    
    
    # make the grids and the profit flow
    grids = makegrids(pea);
    ut = utility(grids,pea);
    tmat = grids.t;
    errcode = 0

    # Initialize all of the convergence criterion variables.  
    in_a_row        = zeros(sz.maxpolit,1); # A vector that checks to see how many times in a row the policy function has converged. 
    sum_in_a_row    = 1000;                 # sum of the last maxpolit policy function differences 
    gap             = 1000.0;               # value function difference
    pgap            = 1000;                 # policy function difference

    # do the VFI
    for iter in 1:sz.maxiter;

        # =====take expectation
        queue = zeros(size(v))

        # Loop over the third dimension 'd'
        for id in 1:sz.nd
            # For each 'd', multiply the transition matrix by the 'np' by 'na' slice of 'v'
            queue[:, :, id] = tmat * v[:, :, id]
        end
 
        # =====interpolate
        if sz.na == sz.npa && sz.nd == sz.npd;
            queuelong = queue;
        else;
            queuelong = fillin(queue,grids);
        end;

        # =====maximize
        # If it's early or if you can do Howard...
        if iter <= sz.earlyiter || sum_in_a_row == 0; 
        #if iter <= 20000 || sum_in_a_row == 0; 

            # If you cannot do Howard 
            if sum_in_a_row > 0; 
                vnew, gidx = maxbellman(queuelong,ut);
            
            # Otherwise do Howard
            else;
                vnew, gidx = howard(queuelong,ut,gidx);
            end;
        # If you are past earlyiter and the policy function has not yet converged.
        else;                           
            vnew, gidx = tinybellman(queuelong,ut,gidx)
        end

        # =====update -- McQueen-Porteus
        vdiff = vnew - v;
        pgap = sum(abs.(gidx.a-old.a)) + sum(abs.(gidx.d-old.d));
        adjfactor  = 0.5*( (beta/(1.0-beta)) * minimum(vdiff) + (beta/(1.0-beta)) * maximum(vdiff));
        gap = maximum(abs.(vdiff));
        in_a_row[2:size(in_a_row,1)] = in_a_row[1:size(in_a_row,1)-1]
        in_a_row[1,1] = pgap
        sum_in_a_row = sum(in_a_row)  

        v = vnew .+ adjfactor;
        old.a = gidx.a; 
        old.d = gidx.d; 
        if mod(iter,10) == 0 && settings.verbose;
            println("iter = ",iter," function diff = ",gap," policy diff = ", pgap)
        end
        if iter == sz.maxiter;
            errcode = 5;
        elseif gap > 10e12;
            errcode = 2;
        end

        if gap < sz.toler || errcode > 0;
            if errcode > 0;
                vnew = 0.0;
                gidx.a = 0;
                gidx.d = 0;
                pol.a .= 0.0;
                pol.d .= 0.0;
                break;
            else; 
                #make the non-integer policy function
                pol.a = makepol(gidx.a,grids.ap);
                pol.d = makepol(gidx.d,grids.dp);
                break;
            end;
        end;
    end;



    outtuple = (v=vnew::Array{Float64}, gidx=gidx, pol=pol, g = grids::NamedTuple, e=errcode::Int);
   
  #  if errcode == 0;
   #     mew = makemew(outtuple); 
   # end;

    #outtuple = (outtuple..., mew=mew::Array{Float64,3})



    return outtuple;
end;
