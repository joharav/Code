module sz;
    const nd        = 11;            #number of points in the durable state grid21
    const na        = 11;            #number of points in the asset state grid21
    const ne        = 4;            #number of points in the exchange rate grid4
    const ny        = 4;            #number of points in the idiosyncratic income grid4
    const npd       = 23;            #number of points in the durable policy grid51
    const npa       = 23;            #number of points in the asset policy grid51
    const pad       = 11;            #number of points to search around the previous point
    const nYears    = 5000;         #number of years to simulate
    const burnin    = 250;           #number of initial years to toss
    const nFirms    = 1000;           #number of firms to simulate
    const nmom      = 6;             #number of moments to calculate
    const maxiter   = 5000;          #maximum number of VFI iterations
    const maxditer  = 1000;          #maximum number of iterations for the stationary distribution
    const distol    = 0.0001;      #tolerance for the stationary distribution
    const toler     = 0.0001;      #VFI tolerance
    const earlyiter = 2;             # Number of times to do the full grid search before you do local searches. 
    const maxpolit  = 5 ;            # Number of times the policy function converges before I believe it.
    const nop       = 17;             #maximum number of parameters
    const noestp    = 5;             #number of parameters you actual estimate 
    const nstd_e    = 3.0;           #number of standard deviations for Tauchen exchange rate, 1.2 for dollar
    const nstd_y    = 1.0;           #number of standard deviations for Tauchen exchange rate, 1.2 for dollar
    const pick = [1,2,3,4,5,6];  # moments to use
end;


module kst;
    const DATA_DIR      = "Data";
    const OUT_DIR       = "Output";
    const MOMS_FILE     = joinpath(DATA_DIR, "datamoments.txt")  ;  # from EFHU
    const W_FILE         = joinpath(DATA_DIR, "W_unclustered.txt") ; # unclustered weighting matrix
    const MNAME_FILE    = joinpath(DATA_DIR, "EFHU_mom_names.txt"); # your moment names
    const PNAME_FILE    = joinpath(DATA_DIR, "EFHU_param_names.txt"); # parameter names
    const PROGRESS_FILE = joinpath(OUT_DIR, "progress.txt");
    const EST_FILE      = joinpath(OUT_DIR, "estfil.txt");
    const RESULTS_FILE  = joinpath(OUT_DIR, "results.txt");

end;

module dtp;

    export Ipol, Pol;


    mutable struct Ipol
        a::Array{Int,5}      # foreign-asset index (a′)
        aa::Array{Int,5}     # local-asset index (aa′)
        d::Array{Int,5}      # durable index (d′)
    end
    
    mutable struct Pol
        aa::Array{Float64,5}  # local-asset policy (levels)
        a::Array{Float64,5}   # foreign-asset policy (levels)
        d::Array{Float64,5}
        c::Array{Float64,5}
    end
end





module globals;
    #  MAKE ABSOLUTELY AS FEW GLOBALS AS POSSIBLE. 
    #set the seed!! 
    using Random, Main.sz;
    Random.seed!(1924);
    draws = rand(sz.nYears+1,sz.nFirms);
end;

module settings; 
    const compstat      = false; 
    const verbose       = false; 
    const irfsshock     = false;     
    const specif_two    = false;
    const welfare       = false;
    const complicated   = false;     # Complicated = true is an identity weight matrix. Otherwise, optimal weight matrix. 

end
