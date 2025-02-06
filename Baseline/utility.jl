function utility(grids::NamedTuple,pea::Vector{Float64}) 

    # Extract grids
    a = grids.a;       # Asset grid
    d = grids.d;       # Durable goods grid
    ap = grids.ap;     # Future asset grid
    dp = grids.dp;     # Future durable goods grid
    pd = grids.p;     # Durable price grid

    # Model parameters
    delta   = pea[2];   # Depreciation rate for durables
    nu      = pea[5];   # Share parameter for nondurables
    gamma   = pea[6];   # Risk aversion
    f       = pea[7];   # Adjustment cost 
    rr      = pea[8];    # Interest rate
    w       = pea[9];    # exogenous income

    # Initialize utility array
    util = zeros(sz.np, sz.na, sz.nd, sz.npa, sz.npd);

    Threads.@threads for iid in 1:sz.npd;
        Threads.@threads for iia in 1:sz.npa;
            Threads.@threads for id in 1:sz.nd;
                Threads.@threads for ia in 1:sz.na;
                    Threads.@threads for ip in 1:sz.np;
                            c = w + a[ia] * (1 + rr) + pd[ip] * (1 - f) * (1 - delta) * d[id] - ap[iia] - pd[ip] * dp[iid]

                            # Check feasibility of consumption and durable goods stock
                            if c > 0 && dp[iid] > 0
                                # Calculate utility
                                util[ip, ia, id, iia, iid] = (((c^nu) * (dp[iid]^(1 - nu)))^(1 - gamma)) / (1 - gamma)
                            else
                                # Apply penalty for infeasible choices
                                util[ip, ia, id, iia, iid] = -1e10
                            end;
            
                    end;
                end;
            end
        end;
    end; 

    return util::Array{Float64};
end

