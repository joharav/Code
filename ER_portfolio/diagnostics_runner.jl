# diagnostics_runner.jl
using DelimitedFiles
using Random, Statistics, Printf, LinearAlgebra, Serialization, Distributions

# 1) Bring in your core code (whatever you usually include before grid search)
include("durable_mod.jl")
include("collectfunctions.jl")          # your big aggregator
include("gmmfunctions_broad.jl")        # objective + moments hooks
include("momentgen.jl")                 # if not already included above
include("diagnostics.jl")               # profiles/slices/sens/sanity
include("ptrue.jl")                     # for sensible defaults

using Main.sz, Main.kst, Main.settings, Main.dtp, Main.globals

# 2) Load θ̂ (your current best). If you have Output/estfil.txt with 5 rows, read it.
function load_theta_hat()
    est_file = joinpath(Main.kst.OUT_DIR, "estfil.txt")
    if isfile(est_file)
        θ = vec(collect(readdlm(est_file)))
        return Float64.(θ)
    else
        # order must match PAR_NAMES in diagnostics.jl
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

    # A) Profiles (dense 1D) — collect() to avoid Dict invariance issues
    grids = Dict{Symbol,AbstractVector}(
        :F_d     => collect(range(0.00, 0.08, length=61)),
        :F_a     => collect(range(0.00, 0.08, length=61)),
        :kappa_a => collect(range(0.20, 1.50, length=61)),
        :chi     => collect(range(0.10, 0.90, length=61)),
        :F_t     => collect(range(0.00, 0.10, length=61))
    )
    grids = Dict{Symbol,AbstractVector}(k => v for (k,v) in grids if haskey(PAR_INDEX, k))
    
    

    run_profiles(θhat; grids=grids, outcsv_dir="Diag/profiles")

    # B) Coarse 2D slices (25×25) for suspect interactions
    pairs = Tuple{Symbol,AbstractVector,Symbol,AbstractVector}[]
    if haskey(PAR_INDEX,:F_d) && haskey(PAR_INDEX,:kappa_a)
        push!(pairs, (:F_d, collect(range(0,0.08,25)),
                       :kappa_a, collect(range(0.2,1.5,25))))
    end
    if haskey(PAR_INDEX,:F_d) && haskey(PAR_INDEX,:F_t)
        push!(pairs, (:F_d, collect(range(0,0.08,25)),
                       :F_t, collect(range(0,0.10,25))))
    end
    run_slices(θhat; pairs=pairs, outcsv_dir="Diag/slices")

    # C) Moment sensitivities (Δm/Δθ bars)
    sens_grids = Dict(
        :F_d     => collect(range(0.00, 0.08, length=41)),
        :kappa_a => collect(range(0.20, 1.50, length=41)),
        :chi     => collect(range(0.10, 0.90, length=41)),
        :F_t     => collect(range(0.00, 0.10, length=41))
    )
    sens_grids = Dict(k=>v for (k,v) in sens_grids if haskey(PAR_INDEX, k))
    run_moment_sensitivities(θhat; grids=sens_grids, outcsv_dir="Diag/sens")

    # D) Sanity check at θ̂
    sc = check_sanity(θhat)
    open("Diag/sanity_θhat.txt", "w") do io
        @printf(io, "ok=%s  adj_rate=%.4f  owner_share=%.4f  usd_asset_share=%.4f\n",
                string(sc.ok), sc.adj_rate, sc.owner_share, sc.usd_asset_share)
    end

    # (Optional) Bounding models quick compare
    # evaluate_bounding_models(θhat)
end
