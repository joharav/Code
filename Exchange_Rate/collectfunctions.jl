# This contains all of the programs to do the VFI homework
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
 include("inbetween.jl")
 include("maxbellman.jl")
 include("howard.jl")
 include("makepol.jl")
 include("mew.jl")
 include("tinybellman.jl")
 include("ptrue.jl")

# ============ VFI =========================================
 include("valfun_adjust.jl")
 include("valfun_noadjust.jl")
 include("valfun.jl")

# ============ Plotting and simulation =====================
 include("plotstuff.jl")
 include("printstuff.jl")
 include("simmodel.jl")
 include("interpol.jl")

# ============Making moments ===============================
 include("momentgen.jl"); #wrapper
 include("makemoments.jl");
