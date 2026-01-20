# ==========================================================================
# Rare-disaster exchange-rate transition (COLUMN-stochastic, grid-preserving)
# Use this with your existing makegrids(ppp) that returns:
#   g.ex :: Vector{Float64}                  (exchange-rate grid, levels)
#   g.te :: Matrix{Float64} size ne×ne       (e-transition, columns sum to 1)
#   g.ty :: Matrix{Float64} size ny×ny       (y-transition, columns sum to 1)
#   g.t  :: Matrix{Float64} size (ne*ny)×(ne*ny) (joint, columns sum to 1)
# ==========================================================================

# ---------- probability mapping ----------
@inline function annual_to_period_prob(pi_y::Float64, periods_per_year::Int)
    @assert 0.0 ≤ pi_y ≤ 1.0
    @assert periods_per_year ≥ 1
    return 1.0 - (1.0 - pi_y)^(1.0 / periods_per_year)
end

# ---------- nearest index on a sorted grid ----------
@inline function nearest_index(grid::AbstractVector{<:Real}, x::Real)
    j = searchsortedlast(grid, x)
    if j < 1
        return 1
    elseif j >= length(grid)
        return length(grid)
    else
        return (abs(x - grid[j]) ≤ abs(grid[j+1] - x)) ? j : (j + 1)
    end
end

# ---------- build map: ie_next -> ie_next_if_jump on LOG grid ----------
function build_jump_map_ex(ex_grid::Vector{Float64}, kappa_e_log::Float64)
    @assert all(ex_grid .> 0.0) "ex_grid must be strictly positive to take logs"
    logex = log.(ex_grid)
    ne = length(ex_grid)
    jm = Vector{Int}(undef, ne)
    @inbounds for iep in 1:ne
        jm[iep] = nearest_index(logex, logex[iep] + kappa_e_log)
    end
    return jm
end

"""
    make_te_disaster_from_base(te_base, ex_grid; pi_annual, periods_per_year, kappa_e_log)

Given baseline e-transition te_base (COLUMN-stochastic: columns sum to 1),
construct disaster mixture:

    te_dis = (1-π)*te_base + π*te_jump

where te_jump remaps next-state mass via a log-jump of size kappa_e_log:
mass in row ie_next is moved to row jump_map[ie_next] (within each column).

Returns:
- te_dis::Matrix{Float64} (column-stochastic)
- meta::NamedTuple(π, jump_map, te_jump)
"""
function make_te_disaster_from_base(te_base::AbstractMatrix{<:Real},
                                    ex_grid::Vector{Float64};
                                    pi_annual::Float64,
                                    periods_per_year::Int,
                                    kappa_e_log::Float64)

    ne1, ne2 = size(te_base)
    @assert ne1 == ne2 "te_base must be square"
    @assert ne1 == length(ex_grid) "te/ex size mismatch"

    te = Float64.(te_base)

    # enforce column-stochastic (your code uses prob = te[ie_next, ie])
    colerr = maximum(abs.(sum(te, dims=1) .- 1.0))
    @assert colerr < 1e-10 "te_base must have columns summing to 1 (max col err = $colerr)"

    π = annual_to_period_prob(pi_annual, periods_per_year)
    jump_map = build_jump_map_ex(ex_grid, kappa_e_log)

    te_jump = zeros(Float64, ne1, ne2)

    # Column-stochastic remap: for each current state (column ie), move mass across rows
    @inbounds for ie in 1:ne2
        for iep in 1:ne1
            p = te[iep, ie]
            p == 0.0 && continue
            iepJ = jump_map[iep]
            te_jump[iepJ, ie] += p
        end
    end

    te_dis = (1.0 - π) .* te .+ π .* te_jump

    # numerical renormalization (columns)
    te_dis ./= sum(te_dis, dims=1)

    return te_dis, (π=π, jump_map=jump_map, te_jump=te_jump)
end

"""
    make_joint_t(te, ty)

Both te and ty must be COLUMN-stochastic.
Returns joint transition for z=(e,y) stacked as iz = (iy-1)*ne + ie (ie fastest).
Then:
    prob = t[iz_next, iz]
"""
function make_joint_t(te::AbstractMatrix{<:Real}, ty::AbstractMatrix{<:Real})
    colerr_e = maximum(abs.(sum(te, dims=1) .- 1.0))
    colerr_y = maximum(abs.(sum(ty, dims=1) .- 1.0))
    @assert colerr_e < 1e-10 "te must be column-stochastic (max col err = $colerr_e)"
    @assert colerr_y < 1e-10 "ty must be column-stochastic (max col err = $colerr_y)"
    return kron(Float64.(ty), Float64.(te))
end

"""
    makegrids_disaster(ppp; pi_annual, kappa_e_log, periods_per_year=4)

Build baseline grids via makegrids(ppp), then replace ONLY:
- g.te (exchange-rate transition)
- g.t  (joint transition)

Keeps the same exchange-rate grid g.ex (levels). This preserves comparability
across counterfactuals while changing tail-risk beliefs.

Assumes baseline makegrids(ppp) returns at least:
- g.ex :: Vector{Float64}
- g.te :: Matrix{Float64} (ne×ne, column-stochastic)
- g.ty :: Matrix{Float64} (ny×ny, column-stochastic)
- g.t  :: Matrix{Float64} ((ne*ny)×(ne*ny), column-stochastic)
"""
function makegrids_disaster(ppp::Vector{Float64};
                            pi_annual::Float64,
                            kappa_e_log::Float64,
                            periods_per_year::Int=4)

    g0 = makegrids(ppp)

    @assert hasproperty(g0, :ex) "baseline grids must contain g.ex"
    @assert hasproperty(g0, :te) "baseline grids must contain g.te (ne×ne, column-stochastic)"
    @assert hasproperty(g0, :ty) "baseline grids must contain g.ty (ny×ny, column-stochastic)"

    ex0 = g0.ex
    te0 = g0.te
    ty0 = g0.ty

    te_dis, meta = make_te_disaster_from_base(te0, ex0;
        pi_annual=pi_annual,
        periods_per_year=periods_per_year,
        kappa_e_log=kappa_e_log
    )

    t_dis = make_joint_t(te_dis, ty0)

    # return same structure as g0, but with te/t replaced
    g_dis = (; g0..., te = te_dis, t = t_dis)

    return g_dis, meta
end
