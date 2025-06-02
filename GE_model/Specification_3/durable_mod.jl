module sz;
    const nd        = 21;            #number of points in the durable state grid
    const na        = 21;            #number of points in the asset state grid
    const ne        = 7;            #number of points in the exchange rate grid
    const npd       = 51;            #number of points in the durable policy grid
    const npa       = 51;            #number of points in the asset policy grid
    const pad       = 11;            #number of points to search around the previous point
    const nYears    = 10100;         #number of years to simulate
    const burnin    = 250;           #number of initial years to toss
    const nFirms    = 1000;           #number of firms to simulate
    const nmom      = 13;             #number of moments to calculate
    const maxiter   = 5000;          #maximum number of VFI iterations
    const maxditer  = 1000;          #maximum number of iterations for the stationary distribution
    const distol    = 0.0001;      #tolerance for the stationary distribution
    const toler     = 0.0001;      #VFI tolerance
    const earlyiter = 2;             # Number of times to do the full grid search before you do local searches. 
    const maxpolit  = 5 ;            # Number of times the policy function converges before I believe it.
    const nop       = 13;             #maximum number of parameters
    const noestp    = 2;             #number of parameters you actual estimate 
    const nstd_e    = 3.0;           #number of standard deviations for Tauchen exchange rate, 1.2 for dollar
    const pick      = [1; 1; 1; 0; 1];     #moments to use
end;


module dtp;

    export Ipol, Pol;

    mutable struct Ipol
        a::Array{Int,3}
        d::Array{Int,3}
    end

    mutable struct Pol
        a::Array{Float64,3}
        d::Array{Float64,3}
        c::Array{Float64,3}
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
    const verbose       = true; 
    const irfsshock     = false;     
    const specif_two    = false;
    const specif_three  = false;
end
