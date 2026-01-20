module sz;
    # ==========================================================================
    # 4D MODEL: States are (e, y, w, d) where w = aa + e*a (total liquid wealth)
    # Dollar share s = e*a / w is chosen within-period (static portfolio choice)
    # ==========================================================================
    
    # State grids - can be finer now with 4D instead of 5D!
    const nd        = 11;           # durable state grid (was 9)
    const nw        = 15;           # total wealth state grid (NEW - replaces na×na)
    const ne        = 7;            # exchange rate grid (was 5)
    const ny        = 7;            # idiosyncratic income grid (was 5)
    
    # Policy grids
    const npd       = 21;           # durable policy grid (was 19)
    const npw       = 21;           # total wealth policy grid (NEW)
    const ns        = 15;           # dollar share grid for within-period choice (NEW)
    
    # Search parameters
    const pad       = 7;            # local search radius
    
    # Simulation parameters
    const nYears    = 1500;         # years to simulate
    const burnin    = 250;          # initial years to discard
    const nFirms    = 500;          # households to simulate (increased from 400)
    
    # Moments
    const nmom      = 6;            # number of moments
    
    # VFI parameters
    const maxiter   = 10000;        # max VFI iterations
    const maxditer  = 1000;         # max distribution iterations
    const distol    = 0.0001;       # distribution tolerance
    const toler     = 0.0001;       # VFI tolerance
    const earlyiter = 3;            # full grid search iterations
    const maxpolit  = 5;            # policy stability requirement
    
    # Parameter counts
    const nop       = 17;           # total parameters
    const noestp    = 5;            # estimated parameters
    
    # Tauchen discretization
    const nstd_e    = 3.0;          # std devs for exchange rate
    const nstd_y    = 1.0;          # std devs for income
    
    # Moments to use
    const pick = [1,2,3,4,5,6];
end;


module kst;
    const DATA_DIR      = "Data";
    const OUT_DIR       = "Output";
    const MOMS_FILE     = joinpath(DATA_DIR, "datamoments.txt");
    const W_FILE        = joinpath(DATA_DIR, "W_unclustered.txt");
    const MNAME_FILE    = joinpath(DATA_DIR, "EFHU_mom_names.txt");
    const PNAME_FILE    = joinpath(DATA_DIR, "EFHU_param_names.txt");
    const PROGRESS_FILE = joinpath(OUT_DIR, "progress.txt");
    const EST_FILE      = joinpath(OUT_DIR, "estfil.txt");
    const RESULTS_FILE  = joinpath(OUT_DIR, "results.txt");
end;


module dtp;
    export Ipol, Pol;
    
    # Policy indices for 4D model: (w', d', s)
    mutable struct Ipol
        w::Array{Int,4}      # total wealth index (w')
        d::Array{Int,4}      # durable index (d')
        s::Array{Int,4}      # dollar share index (s) - within-period choice
    end
    
    # Policy levels for 4D model
    mutable struct Pol
        w::Array{Float64,4}  # total wealth policy (levels)
        d::Array{Float64,4}  # durable policy (levels)
        s::Array{Float64,4}  # dollar share policy (levels, 0 to 1)
        c::Array{Float64,4}  # consumption (implied)
        aa::Array{Float64,4} # local assets (derived: w*(1-s))
        a::Array{Float64,4}  # dollar assets (derived: w*s/e)
    end
end


module globals;
    using Random, Main.sz;
    Random.seed!(1924);
    draws = rand(sz.nYears+1, sz.nFirms);
end;


module settings; 
    const compstat      = false; 
    const verbose       = true; 
    const irfsshock     = false;     
    const specif_two    = false;
    const welfare       = false;
    const complicated   = false;
end
