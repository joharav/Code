function ptrue(enn::Int64)
    pea = zeros(enn);
    pea[1] = 0.95;                               # beta     = Discount factor
    pea[2] = 0.20;                               # delta    = Depreciation rate
    pea[3] = 0.44;                               # rho_e    = persistence
    pea[4] = 0.64;                               # sigma_e  = exchange rate innovation standard deviation
    pea[5] = 0.88;                               # nu       = Share parameter for nondurable consumption 
    pea[6] = 2.00;                               # gamma    = Risk aversion parameter 
    pea[7] = 0.50;                               # f        = adj cost
    pea[8] = 1;                                  # w        = wage rate
    pea[9] = 0.4;                                # chi      = Required maintenance 
    pea[10] = 5;                                 # pd      = price durable goods
    pea[11] = 0.70;                               # ft      = fixed cost on wage rate
    pea[12] = 0.25;                               # tau      = tax rate
    pea[13] = 1/3;                               # h        = hours worked

    return pea::Vector{Float64};
end    