# ==========================================================================
# 4D MODEL: collectfunctions.jl
# Loads the model code AFTER durable_mod.jl has been included exactly once.
# ==========================================================================

using Printf
using Statistics
using LinearAlgebra

# durable_mod.jl must already have run in this session:
# include("durable_mod.jl")
using Main.sz, Main.kst, Main.dtp, Main.globals, Main.settings

# -----------------------------
# Params / grids / tauchen
# -----------------------------
include("ptrue.jl")
include("tauchen.jl")
include("makegrids.jl")

# -----------------------------
# Interpolation
# -----------------------------
include("interpol.jl")

# -----------------------------
# Utility
# -----------------------------
include("utility.jl")
include("utility_noadjust.jl")

# -----------------------------
# Bellman / policy / VFI
# -----------------------------
include("maxbellman.jl")
include("maxbellman_noadjust.jl")
include("tinybellman.jl")
include("tinybellman_noadjust.jl")
include("howard.jl")
include("howard_noadjust.jl")

include("fillin.jl")
include("makepol.jl")
include("makepol_c.jl") 

# IMPORTANT: load the regime solvers first, wrapper last
include("valfun_adjust.jl")
include("valfun_noadjust.jl")
include("valfun.jl")

# -----------------------------
# Simulation / distribution
# -----------------------------
include("simmodel.jl")
include("simmodel_girf.jl")
include("ergodic.jl")

# -----------------------------
# Moments / wrapper
# -----------------------------
include("makemoments.jl")
include("momentgen.jl")

# -----------------------------
# Welfare / diagnostics / plots
# -----------------------------
include("welfare.jl")
include("diagnostics_plots.jl")
include("plotgaps.jl")

# ==========================================================================
# Grid builder (4D always baseline)
# ==========================================================================
grid_builder(pea::Vector{Float64}) = makegrids(pea)

# ==========================================================================
# Smoke test
# ==========================================================================
function test_4D_model(; p::Vector{Float64}=ptrue(sz.nop))
    println("="^60)
    println("Testing 4D Model")
    println("="^60)
    moms = momentgen(p)
    println("="^60)
    return moms
end
