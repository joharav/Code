module sz;
    const nd        = 21;            #number of points in the durable state grid
    const na        = 21;            #number of points in the asset state grid
    const np        = 11;            #number of points in the shock grid
    const ne        = 11;            #number of points in the exchange rate grid
    const npd       = 51;            #number of points in the durable policy grid
    const npa       = 51;            #number of points in the asset policy grid
    const pad       = 11;            #number of points to search around the previous point
    const nYears    = 10100;         #number of years to simulate
    const burnin    = 100;           #number of initial years to toss
    const nFirms    = 100;           #number of firms to simulate
    const nmom      = 14;             #number of moments to calculate
    const maxiter   = 5000;          #maximum number of VFI iterations
    const maxditer  = 1000;          #maximum number of iterations for the stationary distribution
    const distol    = 0.00001;      #tolerance for the stationary distribution
    const toler     = 0.00001;      #VFI tolerance
    const earlyiter = 2;             # Number of times to do the full grid search before you do local searches. 
    const maxpolit  = 5 ;            # Number of times the policy function converges before I believe it.
    const nop       = 11;             #maximum number of parameters
    const noestp    = 2;             #number of parameters you actual estimate 
    const nstd      = 3.0;           #number of standard deviations for Tauchen, 1.2 for dollar
    const nstd_e    = 4.0;           #number of standard deviations for Tauchen exchange rate, 1.2 for dollar
    const pick      = [1; 1; 1; 0; 1];     #moments to use
end;

module kst;
    const beta   = 0.95;       # Discount factor
    const delta  = 0.1;        # Depreciation rate for durable goods
    const rho    = 0.9;        # AR(1) persistence for durable price
    const rho_e  = 0.7;      # AR(1) persistence for exchange rate
    const sigma  = 0.2;        # Volatility of durable price shock
    const sigma_e= 0.3;      # Volatility of exchange rate shock
    const nu     = 0.88;        # Share parameter for nondurable consumption
    const gamma  = 2;          # Risk aversion parameter
    const f      = 0.1;        # Adjustment fixed cost
    const w      = 100.0;        # Wage rate
    const chi    = 0.80;     # Required maintenance 

end;

module dtp;

    export Ipol, Pol;

    mutable struct Ipol
        a::Array{Int,4}
        d::Array{Int,4}
    end

    mutable struct Pol
        a::Array{Float64,4}
        d::Array{Float64,4}
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
    const compstat      = true; 
    const verbose       = true; 
    const complicated   = false;     # Complicated = true is an identity weight matrix. Otherwise, optimal weight matrix. 
end
