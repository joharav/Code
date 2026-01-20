# diagnostics_runner.jl
using DelimitedFiles
using Random, Statistics, Printf, LinearAlgebra, Serialization, Distributions

# 1) Bring in your core code (whatever you usually include before grid search)
include("durable_mod.jl")
include("collectfunctions.jl")          # your big aggregator
include("gmmfunctions_broad.jl")        # objective + moments hooks
include("momentgen.jl")                 # if not already included above
include("diagnostics.jl")               # profiles/slices/sens/sanity (defines Fd_sym, Ft_sym)
include("ptrue.jl")                     # for sensible defaults

using Main.sz, Main.kst, Main.settings, Main.dtp, Main.globals

# -------------------------------
# === Data & weights (cached) ===
# -------------------------------
const DAT = let
    dm     = vec(collect(readdlm(kst.MOMS_FILE)))
    W0     = collect(readdlm(kst.W_FILE))
    mnames = vec(collect(readdlm(kst.MNAME_FILE)))
    pnames = vec(collect(readdlm(kst.PNAME_FILE)))
    (; dm, W0, mnames, pnames)
end

const CH   = sz.pick
const DM   = DAT.dm[CH]
const WRAW = DAT.W0[CH, CH]
const MN   = DAT.mnames[CH]
const PN   = DAT.pnames

const W = let
    Wsym  = Symmetric((WRAW + WRAW')/2)
    ridge = 1e-10
    inv(Wsym + ridge*I)
end

const BIGPEN     = 1e12
const ALL_MN     = vec(collect(readdlm(kst.MNAME_FILE)))  # full names (pre-pick)
const BUDGET_MIN = 180                    # in-Julia budget; leaves cushion
const T0         = time()
const TDEAD      = T0 + 60.0 * BUDGET_MIN
const EPS        = 1e-12
const seed       = 1924
const n_trials   = 2000

# 2) Load θ̂ (your current best). If you have Output/estfil.txt with 5 rows, read it.
function load_theta_hat()
    est_file = joinpath(Main.kst.OUT_DIR, "estfil.txt")
    if isfile(est_file)
        θ = vec(collect(readdlm(est_file)))
        return Float64.(θ)
    else
        # order must match PAR_NAMES in gmmfunctions_broad.jl / diagnostics.jl:
        # [:nu, Symbol("F^d"), :kappa, :chi, Symbol("F^t")]
        return [0.592873, 0.441462, 0.773552, 0.624235, 0.653064]
    end
end

function main()
    # ensure output dirs exist
    mkpath("Diag")
    mkpath("Diag/profiles")
    mkpath("Diag/slices")
    mkpath("Diag/sens")

    θhat = load_theta_hat()

    # =========================================================
    # A) Profiles (dense-ish 1D) — keys must match PAR_INDEX!!
    # =========================================================
    # Fd_sym and Ft_sym are already const globals from diagnostics.jl

    grids = Dict{Symbol,AbstractVector}(
        :nu    => collect(range(0.35, 0.90, length=31)),
        Fd_sym => collect(range(0.001, 0.80, length=31)),
        :kappa => collect(range(0.001, 0.80, length=31)),
        :chi   => collect(range(0.01, 0.90, length=31)),
        Ft_sym => collect(range(0.01, 0.10, length=31))
    )

    # keep only parameters actually in PAR_INDEX
    grids = Dict{Symbol,AbstractVector}(
        k => v for (k,v) in grids if haskey(PAR_INDEX, k)
    )

    run_profiles(θhat; grids=grids, outcsv_dir="Diag/profiles")

    # ==============================================
    # B) Coarse 2D slices for suspect interactions
    # ==============================================
    pairs = Tuple{Symbol,AbstractVector,Symbol,AbstractVector}[]

    if haskey(PAR_INDEX, Fd_sym) && haskey(PAR_INDEX, :kappa)
        push!(pairs, (
            Fd_sym, collect(range(0.0, 0.08, 25)),
            :kappa, collect(range(0.20, 1.50, 25))
        ))
    end

    if haskey(PAR_INDEX, Fd_sym) && haskey(PAR_INDEX, Ft_sym)
        push!(pairs, (
            Fd_sym, collect(range(0.0, 0.08, 25)),
            Ft_sym, collect(range(0.0, 0.10, 25))
        ))
    end

    run_slices(θhat; pairs=pairs, outcsv_dir="Diag/slices")

    # ========================================
    # C) Moment sensitivities (Δm/Δθ bars)
    # ========================================
    sens_grids = Dict{Symbol,AbstractVector}(
        :nu    => collect(range(0.35, 0.90, length=31)),
        Fd_sym => collect(range(0.001, 0.80, length=31)),
        :kappa => collect(range(0.001, 0.80, length=31)),
        :chi   => collect(range(0.01, 0.90, length=31)),
        Ft_sym => collect(range(0.01, 0.10, length=31))
    )

    sens_grids = Dict{Symbol,AbstractVector}(
        k => v for (k,v) in sens_grids if haskey(PAR_INDEX, k)
    )

    run_moment_sensitivities(θhat; grids=sens_grids, outcsv_dir="Diag/sens")

    # ===============================
    # D) Sanity check at θ̂
    # ===============================
    sc = check_sanity(θhat)
    open("Diag/sanity_θhat.txt", "w") do io
        @printf(io, "ok=%s  adj_rate=%.4f  owner_share=%.4f  usd_asset_share=%.4f\n",
                string(sc.ok), sc.adj_rate, sc.owner_share, sc.usd_asset_share)
    end

    # (Optional) bounding models quick compare
    # evaluate_bounding_models(θhat)
end
