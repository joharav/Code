# ============================================================================================================================
#                                                    MASTER PROGRAM FILE
#                  This file includes all model components: Grids, VFI, Simulation, GE, Plotting, Welfare
# ============================================================================================================================

# ======= Grids, Shocks, Utility Functions ========
include("makegrids.jl")
include("tauchen.jl")
include("utility.jl")
include("utility_noadjust.jl")

# ======= Bellman Iteration Components ============
include("fillin.jl")
include("tinybellman.jl")
include("tinybellman_noadjust.jl")
include("maxbellman.jl")
include("maxbellman_noadjust.jl")
include("howard.jl")
include("howard_noadjust.jl")
include("ptrue.jl")

# ======= VFI & Policy Computation ===============
include("makepol.jl")
include("makepol_c.jl")
include("makepol_d_na.jl")

include("valfun_adjust.jl")
include("valfun_noadjust.jl")
include("valfun.jl")

# ======= Simulation Routines =====================
include("simmodel.jl")
include("simmodel_girf.jl")
include("interpol.jl")
include("decision_rules.jl")
include("d_adjust_time_size.jl")
include("aggregate_series.jl")

# ======= Impulse Responses / Shocks ==============
include("mit_shock.jl")
include("girf.jl")

# ======= Plotting Tools ==========================
include("plotstuff.jl")
include("plotgaps.jl")
include("plotdensities.jl")
include("plotgaps_shock.jl")
include("plot_comparison.jl")
include("printstuff.jl")

# ======= GE and Market Clearing ==================
include("GE_mkt_clearing.jl")

# ======= Moment Generation & Welfare =============
include("momentgen.jl")       # wrapper
include("makemoments.jl")
include("adj_gaps_sim.jl")
include("welfare.jl")
