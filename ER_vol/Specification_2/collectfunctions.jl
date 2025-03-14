# This contains all of the programs
# ============================================================================================================================
# ============================================================================================================================
#      FUNCTIONS FOLLOW
# ============================================================================================================================
# ============================================================================================================================

# ============ Stuff that goes into VFI ====================
 include("../makegrids.jl")
 include("../tauchen.jl")
 include("utility.jl")
 include("utility_noadjust.jl")
 include("../fillin.jl")
 #include("inbetween.jl")
 include("../maxbellman.jl")
 include("../maxbellman_noadjust.jl")
 include("../howard.jl")
include("../howard_noadjust.jl")
 include("../makepol.jl")
 include("makepol_c.jl")
 include("../makepol_d_na.jl")
 include("../mew.jl")
 include("../tinybellman.jl")
 include("../tinybellman_noadjust.jl")
 include("../ptrue.jl")

# ============ VFI =========================================
 include("../valfun_adjust.jl")
 include("../valfun_noadjust.jl")
 include("../valfun.jl")

# ============ Plotting and simulation =====================
 include("plotstuff.jl") 
 include("plotgaps.jl") 
 include("plotdensities.jl")   
 include("printstuff.jl") 
 include("../simmodel.jl") 
 include("girf.jl") 
 include("../simmodel_girf.jl")
 include("../interpol.jl")

# ============Making moments ===============================
 include("../momentgen.jl"); #wrapper
 include("../makemoments.jl");
 #include("adjustment_gaps.jl")
 include("../adj_gaps_sim.jl")
 include("../welfare.jl")


