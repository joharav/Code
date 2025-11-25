# diagnostics.jl — identification maps for your SMM

using LinearAlgebra, Statistics, Printf
using DelimitedFiles

# --- Parameter bookkeeping (names must match your 5-D estimation vector) ---

# define the weird ones once
const Fd_sym = Symbol("F^d")
const Ft_sym = Symbol("F^t")

const PAR_NAMES = [:nu, Fd_sym, :kappa, :chi, Ft_sym]   # from gmmfunctions_broad.jl
const PAR_INDEX = Dict(p => i for (i,p) in enumerate(PAR_NAMES))

# --- Objective and moments hooks (rely on your existing methods) ---
# Objective: call the 2-arg fcn, passing +Inf as "best so far"
obj(p::Vector{Float64}) = fcn(p, Inf)

# Simulated moments in the same order as your data/pick
sim_mom(p::Vector{Float64}) = fcn(p)  # 1-arg method returns moments

# --- Profiles ---------------------------------------------------------------
function profile_loss(θhat::Vector{Float64}; j::Symbol, grid)
    jidx = PAR_INDEX[j]
    vals = similar(collect(grid), Float64)
    @inbounds for (k, x) in pairs(grid)
        θ = copy(θhat); θ[jidx] = x
        vals[k] = obj(θ)
    end
    return vals
end

function run_profiles(θhat; grids::Dict{Symbol,AbstractVector}, outcsv_dir="profiles")
    isdir(outcsv_dir) || mkpath(outcsv_dir)
    for (p, grid) in grids
        vals = profile_loss(θhat; j=p, grid)
        fn = joinpath(outcsv_dir, "$(p)_profile.csv")
        open(fn, "w") do io
            @printf(io, "param,theta,loss\n")
            for (x,v) in zip(grid, vals)
                @printf(io, "%s,%.10f,%.10f\n", String(p), x, v)
            end
        end
    end
end

# --- 2D slices --------------------------------------------------------------
function slice2D(θhat; p1::Symbol, g1, p2::Symbol, g2)
    j1, j2 = PAR_INDEX[p1], PAR_INDEX[p2]
    Z = Matrix{Float64}(undef, length(g1), length(g2))
    for (i, x1) in enumerate(g1), (k, x2) in enumerate(g2)
        θ = copy(θhat); θ[j1] = x1; θ[j2] = x2
        Z[i,k] = obj(θ)
    end
    return Z
end

function run_slices(
    θhat;
    pairs::Vector{Tuple{Symbol,AbstractVector,Symbol,AbstractVector}},
    outcsv_dir="slices"
)
    isdir(outcsv_dir) || mkpath(outcsv_dir)
    for (p1, g1, p2, g2) in pairs
        Z = slice2D(θhat; p1, g1, p2, g2)
        fn = joinpath(outcsv_dir, "$(p1)_$(p2)_slice.csv")
        open(fn, "w") do io
            @printf(io, "p1,p2,loss\n")
            for (i,x1) in enumerate(g1), (k,x2) in enumerate(g2)
                @printf(io, "%.10f,%.10f,%.10f\n", x1, x2, Z[i,k])
            end
        end
    end
end

# --- Moment sensitivities Δm/Δθ -------------------------------------------
function moment_sensitivity(θhat; j::Symbol, grid::AbstractVector)
    jidx = PAR_INDEX[j]
    X = collect(grid)
    Ms = Vector{Vector{Float64}}(undef, length(X))
    for (k, x) in enumerate(X)
        θ = copy(θhat); θ[jidx] = x
        Ms[k] = sim_mom(θ)  # vector of moments (your pick order)
    end
    Mgrid = reduce(vcat, (permutedims(m) for m in Ms))  # |grid| × n_m
    # local slope via OLS on small window around θhat[j]
    at = X[argmin(abs.(X .- θhat[jidx]))]
    idx = sortperm(abs.(X .- at))[1:min(5, length(X))]
    slopes = Vector{Float64}(undef, size(Mgrid,2))
    for m in axes(Mgrid,2)
        xw = X[idx]; yw = Mgrid[idx, m]
        β = [ones(length(idx)) xw] \ yw
        slopes[m] = β[2]
    end
    return (; grid=X, Mgrid, slopes)
end

function run_moment_sensitivities(θhat; grids::Dict{Symbol,AbstractVector}, outcsv_dir="sens")
    isdir(outcsv_dir) || mkpath(outcsv_dir)
    # names from your files (keeps alignment with datamoments pick)
    mnames = vec(collect(readdlm(joinpath(Main.kst.DATA_DIR, "EFHU_mom_names.txt"))))
    pick   = Main.sz.pick
    mnames = mnames[pick]
    for (p, grid) in grids
        res = moment_sensitivity(θhat; j=p, grid)
        # full grid of moments
        fn = joinpath(outcsv_dir, "$(p)_moments_over_grid.csv")
        open(fn, "w") do io
            @printf(io, "theta")
            for nm in mnames
                @printf(io, ",%s", String(nm))
            end
            @printf(io, "\n")
            for i in eachindex(res.grid)
                @printf(io, "%.10f", res.grid[i])
                for m in axes(res.Mgrid,2)
                    @printf(io, ",%.10f", res.Mgrid[i,m])
                end
                @printf(io, "\n")
            end
        end
        # local slopes
        fn2 = joinpath(outcsv_dir, "$(p)_moment_slopes.csv")
        open(fn2, "w") do io
            @printf(io, "moment,dm_d%s\n", String(p))
            for (nm, s) in zip(mnames, res.slopes)
                @printf(io, "%s,%.10f\n", String(nm), s)
            end
        end
    end
end

# --- Sanity checks: inaction band + cross-section shares -------------------
function check_sanity(θ::Vector{Float64}; target_adj_rate=(0.05,0.40),
                      target_owner_share=(0.40,0.90))
    answ = valfun(buildparam(θ))
    simd = simmodel(answ)
    r0 = (Main.sz.burnin-2):Main.sz.nYears
    adj = vec(simd.adjust_indicator[r0, :] .> 0)
    d   = simd.d[r0, :]
    ex  = simd.ex[r0, :]
    a   = simd.a[r0, :]
    aa  = simd.aa[r0, :]

    adj_rate    = mean(adj)
    owner_share = mean(vec(d) .> 0.0)
    a_loc = aa
    a_fx  = ex .* a
    a_eff = a_loc .+ a_fx .+ 1e-12
    usd_asset_share = mean(vec(a_fx ./ a_eff) .> 0.5)  # coarse, but stable

    ok_adj   = target_adj_rate[1] ≤ adj_rate ≤ target_adj_rate[2]
    ok_owner = target_owner_share[1] ≤ owner_share ≤ target_owner_share[2]
    return (ok = ok_adj && ok_owner,
            adj_rate=adj_rate, owner_share=owner_share, usd_asset_share=usd_asset_share)
end
