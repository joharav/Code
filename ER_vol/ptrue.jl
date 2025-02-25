function ptrue(enn::Int64)
    pea = zeros(enn);
    pea[1] = 0.95;                               # beta     = Discount factor
    pea[2] = 0.10;                               # delta    = Depreciation rate
    pea[3] = 0.70;                               # rho_e    = persistence
    pea[4] = 0.30;                               # sigma_e  = exchange rate innovation standard deviation
    pea[5] = 0.88;                               # nu       = Share parameter for nondurable consumption 
    pea[6] = 2.00;                               # gamma    = Risk aversion parameter 
    pea[7] = 0.80;                               # f        = adj cost
    pea[8] = 1;                                # w        = wage rate
    pea[9] = 0.8;                                # chi      = Required maintenance 
    pea[10] = 4;                               # pd      = price durable goods

    return pea::Vector{Float64};
end