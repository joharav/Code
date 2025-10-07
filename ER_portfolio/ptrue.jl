function ptrue(enn::Int64)
    pea = zeros(enn);
    pea[1] = 0.95;                               # beta     = Discount factor
    pea[2] = 0.05;                               # delta    = Depreciation rate
    pea[3] = 0.66;                               # rho_e    = persistence
    pea[4] = 0.25;                               # sigma_e  = exchange rate innovation standard deviation
    pea[5] = 0.51;                               # nu       = Share parameter for nondurable consumption 
    pea[6] = 2.00;                               # gamma    = Risk aversion parameter 
    pea[7] = 0.03;                               # f        = adj cost 0.3 gives you a 3% adj
    pea[8] = 1;                                  # w        = wage rate
    pea[9] = 0.03;                                # r_foreign      = Foreign rate
    pea[10] = 5;                                 # pd      = price durable goods
    pea[11] = 0.5;                              # kappa    = fixed cost of asset adjustment
    pea[12] = 0.0;                              # tau      = tax rate
    pea[13] = 1;                               # h        = hours worked
    pea[14] = 0.9;                              # rho_y    = AR(1) persistence for idiosyncratic income
    pea[15] = 0.2;                               # sigma_y  = Volatility of idiosyncratic income shock
    pea[16] = 1.0;                               # theta   = Share dollar savings


    return pea::Vector{Float64};
end    