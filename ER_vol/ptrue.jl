function ptrue(enn::Int64)
    pea = zeros(enn);
    pea[1] = 0.95;                               # beta     = Discount factor
    pea[2] = 0.10;                               # delta    = Depreciation rate
    pea[3] = 0.90;                               # rho      = persistence 
    pea[4] = 0.70;                               # rho_e    = persistence
    pea[5] = 0.20;                               # sigma    = shock innovation standard deviation
    pea[6] = 0.30;                               # sigma_e  = exchange rate innovation standard deviation
    pea[7] = 0.88;                               # nu       = Share parameter for nondurable consumption 
    pea[8] = 2.00;                               # gamma    = Risk aversion parameter 
    pea[9] = 0.10;                               # f        = adj cost
    pea[10] = 100;                               # w        = wage rate
    pea[11] = 0.8;                               # chi      = Required maintenance 

    return pea::Vector{Float64};
end