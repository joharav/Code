# ==========================================================================
# 4D MODEL: Default/true parameter values
# ==========================================================================

function ptrue(enn::Int64)
    pea = zeros(enn)
    
    # Calibrated parameters
    pea[1] = 0.98          # beta (discount factor)
    pea[2] = 0.05         # delta (depreciation rate)
    pea[3] = 0.66          # rho_e (ER persistence)
    pea[4] = 0.15          # sigma_e (ER volatility)
    pea[6] = 2.0           # gamma (risk aversion)
    pea[8] = 1.0           # wage
    pea[9] = 0.0075          # r_foreign (dollar interest rate)
    pea[10] = 1.0          # pd (durable price)
    pea[12] = 0.0         # tau (tax rate)
    pea[13] = 1.0          # h (labor supply)
    pea[14] = 0.9          # rho_y (income persistence)
    pea[15] = 0.20         # sigma_y (income volatility)
    
    # Estimated parameters - UPDATED for 4D model
    # These are starting points that should match data moments better
    pea[5] = 0.58          # nu (non-durable share) - targets dwealth_mean
    pea[7] = 0.4           # F_d (durable fixed cost) - targets adj_rate, duration
    pea[11] = 0.20         # kappa (dollar transaction cost) - targets dollar_share
    pea[16] = 0.50         # chi (maintenance effectiveness) - targets dwealth_var
    pea[17] = 0.30         # F_t (time cost)
    
    return pea
end


# ==========================================================================
# Parameter bounds for estimation
# ==========================================================================
function param_bounds()
    lb = zeros(sz.noestp)
    ub = zeros(sz.noestp)
    
    # Order: [nu, F_d, kappa, chi, F_t]
    lb[1] = 0.40;  ub[1] = 0.70   # nu
    lb[2] = 0.50;  ub[2] = 5.00   # F_d - MUCH higher to reduce adj rate
    lb[3] = 0.01;  ub[3] = 0.50   # kappa - LOWER to increase dollar share
    lb[4] = 0.30;  ub[4] = 0.70   # chi
    lb[5] = 0.10;  ub[5] = 0.60   # F_t
    
    return lb, ub
end


# ==========================================================================
# Parameter names
# ==========================================================================
const PARAM_NAMES = [
    "beta",        # 1
    "delta",       # 2
    "rho_e",       # 3
    "sigma_e",     # 4
    "nu",          # 5
    "gamma",       # 6
    "F_d",         # 7
    "wage",        # 8
    "r_foreign",   # 9
    "p_d",         # 10
    "kappa",       # 11
    "tau",         # 12
    "h",           # 13
    "rho_y",       # 14
    "sigma_y",     # 15
    "chi",         # 16
    "F_t"          # 17
]

const PARAM_LABELS = Dict(
    "beta"      => "β",
    "delta"     => "δ",
    "rho_e"     => "ρₑ",
    "sigma_e"   => "σₑ",
    "nu"        => "ν",
    "gamma"     => "γ",
    "F_d"       => "Fᵈ",
    "wage"      => "w",
    "r_foreign" => "r*",
    "p_d"       => "pᵈ",
    "kappa"     => "κ",
    "tau"       => "τ",
    "h"         => "h",
    "rho_y"     => "ρᵧ",
    "sigma_y"   => "σᵧ",
    "chi"       => "χ",
    "F_t"       => "Fᵗ"
)
