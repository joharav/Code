function simann(tuning_parameters,x_start,enn,lb,ub,iseed,tee,vee_em);

    # Date: January 27, 2022
    # Author: Toni Whited
    # Translated from the Alan Miller/Bill Goffe Fortran 90 version 
    # Contains many fortranesque contructions that are not pretty in Julia. 

    # Input: 
    #   x_start = the starting value (this is actually a tuple that contains other model constants)
    #   enn     = length of x_start 
    #   tuning_parameters are described below 
    #   lb and ub are upper and lower bounds for the parameters 
    #   iseed is the random number generator seed 
    
    #  Input-Output
    #tee     = answ.t;       The final temperature.
    #vee_em  = answ.v;       The final step size  

    # Output: 
    #xopt    = answ.x;       The parameters that optimize the function. 
    #fopt    = answ.f;       The optimized function value 
    #nacc    = answ.a;       The number of accepted function evaluations.                              
    #nobds   = answ.b;       The total number of trial function evaluations that         
    #                        would have been out of bounds of LB and UB. Note that        
    #                        a trial point is randomly selected between LB and UB.  
    #ier     = answ.i;       The error return number.                         
    #                        Values: 0 - Normal return; termination criteria achieved.  
    #                                1 - Number of function evaluations (NFCNEV) is     
    #                                    greater than the maximum number (MAXEVL).      
    #                                2 - The starting value (X) is not inside the       
    #                                    bounds (LB and UB).                            
    #                                3 - The initial temperature is not positive.       
    #                                99 - Should not be seen; only used internally.     

    # =============Extract all of the tuning parameters ================================
    mini    = tuning_parameters.mini;           #are you maximizing or minimizing 
    rt      = tuning_parameters.rt;             #temperature reduction factor. I use .85 
    eps     = tuning_parameters.eps;            #convergence tolerance. If the final function
                                                #values from the last neps temperatures differ from the
                                                #corresponding value at the current temperature by less than
                                                #EPS and the final function value at the current temperature
                                                #differs from the current optimal function value by less than
                                                #EPS, execution terminates and IER = 0 is returned. (EP)       
    ns      = tuning_parameters.ns;             #number of cycles 
                                                #After NS*N function evaluations, each element of VM is adjusted so 
                                                #that approximately half of all function evaluations
                                                #are accepted.  The suggested value is 20.  
                                                #The higher NS the slower the step size is adjusted.    
    nt      = tuning_parameters.nt;             #Number of iterations before temperature reduction. After
                                                #NT*NS*N function evaluations, temperature (T) is changed
                                                #by the factor RT.  Value suggested by Corana et al. is
                                                #MAX(100, 5*N).  See Goffe et al. for further advice.                        
    neps    = tuning_parameters.neps;           #NEPS - Number of final function values used to decide upon termination.  See EPS.  
                                                #Suggested value is 4. 
    maxevl  = tuning_parameters.maxevl;         #The maximum number of function evaluations.  If it is exceeded, IER = 1. 
    cee     = tuning_parameters.cee;            #Vector that controls the step length adjustment.  The suggested value for all elements is 2.0. 
    iprint  = tuning_parameters.iprint          #Controls how much gets printed out to the output file. 
                                                #1 is just the basics
                                                #2 is a little bit more
                                                #3 is a lot 
    # ============================================End of preamble ===========================================
    
    # ============================================Start of initializations=========================================================================

    # == initialize outputs so that the program will have something to return if it barfs

    ier = 99
    answ = (x=0,f=0,a=0,b=0,i=ier,t=tee,v=vee_em);  

    # =====Stop the program if something is deeply wrong
    if tee < 0;
        ier = 3; 
        println(" The initial temp is negative. Don't do that.")
        return answ; 
    end;
    out_of_bounds = sum(collect(lb .> x_start .|| ub .< x_start));
    if out_of_bounds > 0;
        ier = 2; 
        println(" Starting value out of bounds. Don't do that.")
        return answ; 
    end; 

    # == Open the output file. == 

    simannout = "Output/simannout.txt";
    io = open(simannout,"w")
    println(io,"Simulated Annealing Output")
    close(io)
    
    # == initialize the internal variables that are not scalars
    fstar   = ones(neps).*1.0e20;        #this stores the last neps optimal function values for the last neps temps 
    xp      = zeros(enn);                #candidate parameter value. This gets reinitialized every trial. 
    nacp    = zeros(enn);                #This is the number of acceptances per parameter per cycle
                                         #Reset every cycle. Used to adjust step size.  

    # Set the initial values
    nacc = 0        #number accepted       
    nobds = 0       #number out of bounds                                 
    x = x_start;    #x is the internal parameter value
    xopt = x;
    fopt = 0.0; 

    #set the seed!! 
    Random.seed!(iseed);

    # Initialize the optimal function value 

    f = fcn(x, fopt)   
    if mini;
        f = -f;
    end; 
    fopt = f
    fstar[1] = f  

    # ============================================End of initializations===========================================================================   

    # Start the SA iterations  
    nrej = 0        #number rejected
    nnew = 0        #number new optima
    ndown = 0       #number of down steps
    nup = 0         #number of up steps
    nfcnev = 1      #SA trial
    
    while nfcnev <= maxevl;
        
        println("I am on evaluation number ",nfcnev)

        for m in 1:nt
            
            for j in 1:ns
                
                for h in 1:enn

                    #  If too many function evaluations occur, terminate the algorithm.
                    if (nfcnev >= maxevl)
                        ier = 1
                        answ = (x=xopt,f=fopt,a=nacc,b=nobds,i=ier,t=tee,v=vee_em);
                        return answ;
                    end;
          
                    #  Generate XP, the trial value of X. Note use of VM to choose XP.
                    xp = zeros(enn)
                    for i in 1:enn
                        if (i == h)
                            xp[i] = x[i] + (rand(1)[1]*2. - 1.) * vee_em[i]
                        else;
                            xp[i] = x[i]
                        end;
                  
                        #  If XP is out of bounds, select a point in bounds for the trial.
                        if ((xp[i] < lb[i]) || (xp[i] > ub[i]))
                            xp[i] = lb[i] + (ub[i] - lb[i])*rand(1)[1]
                            nobds = nobds + 1
                        end;
                    end;
          
                    #  Evaluate the function with the trial point XP and return as FP.min 
                    fp = fcn(xp, (-2*Int(mini)+1.)*fopt);   
                    if mini;
                        fp = -fp;
                    end; 
                    nfcnev = nfcnev + 1
          
                    #  Accept the new point if the function value increases.
                    if (fp >= f)
                        
                        #if iprint == 3;
                        #    write(simannout,"  POINT ACCEPTED")
                        #end 
                        x = xp
                        f = fp; 
                        nacc = nacc + 1
                        nacp[h] = nacp[h] + 1
                        nup = nup + 1
                    
                        #  If greater than any other point, record as new optimum.
                        if (fp > fopt)
                            if (iprint == 3)
                                prt8(simannout, vee_em, xopt, x, tee)
                            end
                            xopt = xp
                            fopt = fp;  
                            nnew = nnew + 1
                        end;
          
                    #  If the point is lower, use the Metropolis criteria to decide on
                    #  acceptance or rejection.
                    else;
                        
                        p = exp((fp - f)/tee)    
                        pp = rand(1)[1]
                        if (pp < p)
                            x = xp
                            f = fp ;   
                            nacc = nacc + 1
                            nacp[h] = nacp[h] + 1
                            ndown = ndown + 1
                            #println("went down")
                        else;
                            nrej = nrej + 1
                        end;
                    end;
          
                end; #end of the parameter loop
            end; #end of the cycle loop 
          
            #  Adjust VM so that approximately half of all evaluations are accepted.
            for i in 1:enn; 
                ratio = float(nacp[i]) /float(ns)
                if (ratio > .6)
                    vee_em[i] = vee_em[i]*(1. + cee*(ratio - .6)/.4)
                elseif (ratio < .4)
                    vee_em[i] = vee_em[i]/(1. + cee*((.4 - ratio)/.4))
                end;

                if (vee_em[i] > (ub[i]-lb[i]))
                    vee_em[i] = ub[i] - lb[i]
                end;
                
            end;
          
            if (iprint >= 2)
                prt8(simannout, vee_em, xopt, x, tee)
            end;
          
            nacp = zeros(enn);
          
        end; #end of temperature reduction loop 

        if (iprint >= 1) 
            prt9(simannout, tee, xopt, vee_em, fopt, nup, ndown, nrej, nobds, nnew, nfcnev)
        end;

        #!  Check termination criteria.
        quit_sa = false
        fstar[1] = f
        if ((fopt - fstar[1]) <= eps)
           quit_sa = true
        end

        for ii in 1:neps
            if (abs(f - fstar[ii]) > eps)
                quit_sa = false
            end
        end
      
        #Terminate SA if appropriate.
        if quit_sa
            x = xopt
            f = fopt
            ier = 0
            #if (iprint >= 1)  
            #    write(simannout,"Did it, converged!")
            #end 
            answ = (x=xopt,f=-fopt,a=nacc,b=nobds,i=ier,t=tee,v=vee_em);
            return answ; 
        else
            #If termination criteria is not met, prepare for another loop.
            tee = rt*tee
            for i in neps:-1:2
                fstar[i] = fstar[i-1]
            end;
            f = fopt
            x = xopt
            nnew = 0
        end

    end; #end of final SA convergence loop 

end;

function prt8(simannout, vee_em, xopt, x, t)
    io = open(simannout,"a")
    write(io,"Intermediate results after step length adjustment")
    write(io," ")
    write(io,"new step length (vm)")
    println(io,vee_em)
    write(io," ")
    write(io,"current optimal x")
    println(io,xopt)
    write(io," ")
    write(io,"current  x")
    println(io,x)
    write(io," ")
    write(io,"temperature")
    println(io,t)
    close(io)
   
   filname = "Output/veeem.txt"
   io = open(filname,"w")
   for jj in eachindex(vee_em); @printf(io,"%25.16f\n",vee_em[jj]); end;  
   close(io)
   
   
   filname = "Output/tempfil.txt"
   io = open(filname,"w")
   @printf(io,"%25.16f\n",t);
   close(io)
   end;

   function prt9(simannout, t, xopt, vee_em, fopt, nup, ndown, nrej, nobds, nnew, nfcnev)

    totmov = nup + ndown + nrej
    
    filname = "Output/veeem.txt"
    io = open(filname,"w")
    for jj in eachindex(vee_em); @printf(io,"%25.16f\n",vee_em[jj]); end;
    close(io)
    
    filname = "Output/estfil.txt"
    io = open(filname,"w")
    for jj in eachindex(xopt); @printf(io,"%25.16f\n",xopt[jj]); end;
    close(io)
    
    
    
    filname = "Output/tempfil.txt"
    io = open(filname,"w")
    @printf(io,"%25.16f\n",t);
    close(io)
    
    io = open(simannout,"a")
    @printf(io," intermediate results before next temperature reduction \n")
    @printf(io,"  current temperature:          %25.16f\n", t)
    @printf(io,"  max function value so far:    %25.16f\n", fopt)
    @printf(io,"  total moves:                  %8i\n", totmov  )
    @printf(io,"     uphill:                    %8i\n", nup     )
    @printf(io,"     accepted downhill:         %8i\n", ndown   )
    @printf(io,"     rejected downhill:         %8i\n", nrej    )
    @printf(io,"  out of bounds trials:         %8i\n", nobds  )
    @printf(io,"  new maxima this temperature:  %8i\n", nnew    )
    @printf(io,"  total number of evaluations:  %8i\n", nfcnev  )  
    @printf(io," \n") 
    @printf(io,"new step length (vm) \n")  
    for jj in eachindex(vee_em); @printf(io,"%10.5f",vee_em[jj]); end;
    @printf(io," \n")
    @printf(io," \n")
    @printf(io,"optimal point \n")
    for jj in eachindex(xopt); @printf(io,"%10.5f",xopt[jj]); end;
    @printf(io," \n")   
    @printf(io," \n")   
    @printf(io," \n")   
    close(io)  
end
    
   