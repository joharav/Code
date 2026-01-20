# ==========================================================================
# 4D MODEL: Collect all functions
# Include this file to load the entire model
# ==========================================================================

using Printf
using Statistics
using LinearAlgebra

# Module definitions (must come first)
include("durable_mod.jl")

# Import modules into scope
using .sz
using .kst
using .dtp
using .globals
using .settings

# Parameters
include("ptrue.jl")

# Grid construction
include("tauchen.jl")       # Tauchen discretization (unchanged)
include("makegrids.jl")     # 4D grid construction

# Interpolation
include("interpol.jl")

# Utility functions
include("utility.jl")           # Adjustment regime utility
include("utility_noadjust.jl")  # Non-adjustment utility

# Bellman operators
include("maxbellman.jl")          # Full grid search (adjust)
include("maxbellman_noadjust.jl") # Full grid search (no-adjust)
include("tinybellman.jl")         # Local search (adjust)
include("tinybellman_noadjust.jl") # Local search (no-adjust)
include("howard.jl")              # Howard acceleration (adjust)
include("howard_noadjust.jl")     # Howard acceleration (no-adjust)

# Grid fill-in
include("fillin.jl")

# Policy extraction
include("makepol.jl")

# Value function iteration
include("valfun_adjust.jl")
include("valfun_noadjust.jl")
include("valfun.jl")        # Combined wrapper

# Simulation
include("simmodel.jl")
include("simmodel_girf.jl")  # Impulse response simulation

# Ergodic distribution
include("ergodic.jl")

# Moment computation
include("makemoments.jl")
include("momentgen.jl")

# Welfare
include("welfare.jl")

# Adjustment gaps analysis
include("adj_gaps_sim.jl")

# Plotting utilities (unchanged)
include("plotgaps.jl")

# GMM/SMM functions - include separately in estimation scripts
# include("gmmfunctions.jl")

# ==========================================================================
# Grid builder function (chooses baseline vs disaster grids)
# For 4D model, always use baseline grids
# ==========================================================================
function grid_builder(pea::Vector{Float64})
    return makegrids(pea)
end


# ==========================================================================
# Quick test function
# ==========================================================================
function test_4D_model()
    println("=" ^ 60)
    println("Testing 4D Model")
    println("=" ^ 60)
    
    # Default parameters
    p = ptrue(sz.nop)
    
    println("\nParameters:")
    pnames = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "F_d", "wage",
              "r_star", "pd", "kappa", "tau", "h", "rho_y", "sigma_y", "chi", "F_t"]
    for (i, (name, val)) in enumerate(zip(pnames, p))
        @printf("  %2d. %-10s = %.4f\n", i, name, val)
    end
    
    println("\nRunning model...")
    moms = momentgen(p)
    
    println("\n" * "=" ^ 60)
    println("Test complete!")
    println("=" ^ 60)
    
    return moms
end
# Grid builder function (chooses baseline vs disaster grids)
# For 4D model, always use baseline grids
# ==========================================================================
function grid_builder(pea::Vector{Float64})
    return makegrids(pea)
end


# ==========================================================================
# Quick test function
# ==========================================================================
function test_4D_model()
    println("=" ^ 60)
    println("Testing 4D Model")
    println("=" ^ 60)
    
    # Default parameters (from your estimation)
    p = zeros(17)
    p[1] = 0.98      # beta
    p[2] = 0.025     # delta (depreciation)
    p[3] = 0.9       # rho_e (ER persistence)
    p[4] = 0.15      # sigma_e (ER volatility)
    p[5] = 0.6       # nu (non-durable share)
    p[6] = 2.0       # gamma (risk aversion)
    p[7] = 0.5       # f (durable fixed cost) - INCREASED
    p[8] = 1.0       # wage
    p[9] = 0.02      # r_star (dollar rate)
    p[10] = 1.0      # pd (durable price)
    p[11] = 0.1      # kappa (dollar transaction cost) - DECREASED
    p[12] = 0.25     # tau (tax rate)
    p[13] = 1.0      # h (labor supply)
    p[14] = 0.9      # rho_y (income persistence)
    p[15] = 0.1      # sigma_y (income volatility)
    p[16] = 0.5      # chi (maintenance effectiveness)
    p[17] = 0.1      # ft (time cost)
    
    println("\nParameters:")
    pnames = ["beta", "delta", "rho_e", "sigma_e", "nu", "gamma", "f", "wage",
              "r_star", "pd", "kappa", "tau", "h", "rho_y", "sigma_y", "chi", "ft"]
    for (i, (name, val)) in enumerate(zip(pnames, p))
        @printf("  %2d. %-10s = %.4f\n", i, name, val)
    end
    
    println("\nRunning model...")
    moms = momentgen(p)
    
    println("\n" * "=" ^ 60)
    println("Test complete!")
    println("=" ^ 60)
    
    return moms
end
