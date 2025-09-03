# This contains all of the programs
# ============================================================================================================================
# ============================================================================================================================
#      FUNCTIONS FOLLOW
# ============================================================================================================================
# ============================================================================================================================

# ============ Stuff that goes into VFI ====================
include("makegrids.jl")
include("tauchen.jl")
include("utility.jl")
include("utility_noadjust.jl")
include("fillin.jl")
include("maxbellman.jl")
include("maxbellman_noadjust.jl")
include("howard.jl")
include("howard_noadjust.jl")
include("makepol.jl")
include("makepol_c_twoasset.jl")
include("makepol_d_na.jl")
include("tinybellman.jl")
include("tinybellman_noadjust.jl")
include("ptrue.jl")

# ============ VFI =========================================
include("valfun_adjust.jl")
include("valfun_noadjust.jl")
include("valfun.jl")

# ============ Plotting and simulation =====================
include("plotstuff.jl")
include("plotgaps.jl")
include("plotdensities.jl")
include("printstuff.jl")
include("simmodel.jl") 
include("girf.jl")
include("simmodel_girf.jl")
include("interpol.jl")
include("aggregate_series.jl")
include("d_adjust_time_size.jl") 
include("decision_rules.jl") 
include("plotgaps_shock.jl")
include("plot_ergodic.jl")
include("plot_distribution_panels.jl")
# ============ GMM and SMM =================================
include("gmmfunctions_broad.jl")
#include("inflnc_functions.jl")
#include("inflnc.jl")
include("simann.jl")
#include("empirical_dist.jl") # Not used in the current setup

# ============Making moments ===============================
include("momentgen.jl"); #wrapper
include("makemoments_smm.jl");
include("adj_gaps_sim.jl")
include("welfare.jl")
include("welfare_dispersion.jl")
include("ergodic.jl")
include("welfare_cases.jl")



