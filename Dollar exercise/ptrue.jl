function ptrue(enn::Int64)
    pea = zeros(enn);
    pea[1] = 0.95;                               # beta     = Discount factor
    pea[2] = 0.10;                               # delta    = Depreciation rate
    pea[3] = 0.90;                               # rho      = persistence 
    pea[4] = 0.90;                               # rho_e    = persistence
    pea[5] = 0.20;                               # sigma    = shock innovation standard deviation
    pea[6] = 0.40;                               # sigma_e    = exchange rate innovation standard deviation
    pea[7] = 0.50;                               # nu       = Share parameter for nondurable consumption 
    pea[8] = 2.00;                               # gamma    = Risk aversion parameter 
    pea[9] = 0.20;                               # adj cost
  
    return pea::Vector{Float64};
end