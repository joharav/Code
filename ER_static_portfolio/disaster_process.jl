# ------------------------------------------------------------
# Rare-disaster ER transition (ROW-stochastic)
# te[i, j] = Pr(e_{t+1}=j | e_t=i), rows sum to 1
# ------------------------------------------------------------

@inline function annual_to_period_prob(pi_y::Float64, periods_per_year::Int)
    @assert 0.0 ≤ pi_y ≤ 1.0
    return 1.0 - (1.0 - pi_y)^(1.0 / periods_per_year)
end

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

function build_jump_map_ex(ex_grid::Vector{Float64}, kappa_e_log::Float64)
    @assert all(ex_grid .> 0.0)
    logex = log.(ex_grid)
    ne = length(ex_grid)
    jm = Vector{Int}(undef, ne)
    @inbounds for j in 1:ne
        jm[j] = nearest_index(logex, logex[j] + kappa_e_log)
    end
    return jm
end

"""
    make_te_disaster_from_base_row(te_base, ex_grid; pi_annual, periods_per_year, kappa_e_log)

Inputs/outputs are ROW-stochastic.
Disaster: with prob π, next-state mass is remapped by a log jump of size kappa_e_log.
"""
function make_te_disaster_from_base_row(te_base::AbstractMatrix{<:Real},
                                        ex_grid::Vector{Float64};
                                        pi_annual::Float64,
                                        periods_per_year::Int,
                                        kappa_e_log::Float64)

    ne1, ne2 = size(te_base)
    @assert ne1 == ne2 == length(ex_grid)

    te = Float64.(te_base)

    # check row-stochastic
    rowerr = maximum(abs.(sum(te, dims=2) .- 1.0))
    @assert rowerr < 1e-10 "te_base must be ROW-stochastic (max row err=$rowerr)"

    π = annual_to_period_prob(pi_annual, periods_per_year)
    jump_map = build_jump_map_ex(ex_grid, kappa_e_log)

    te_jump = zeros(Float64, ne1, ne2)

    # For each current state i (row), move probability mass across NEXT states (columns)
    @inbounds for i in 1:ne1
        for j in 1:ne2
            p = te[i, j]
            p == 0.0 && continue
            jJ = jump_map[j]
            te_jump[i, jJ] += p
        end
    end

    te_dis = (1.0 - π) .* te .+ π .* te_jump

    # renormalize rows
    te_dis ./= sum(te_dis, dims=2)

    return te_dis, (π=π, jump_map=jump_map, te_jump=te_jump)
end

"""
    make_joint_t_row(te, ty)

Both ROW-stochastic. Returns ROW-stochastic joint on (e,y) with e fastest:
iz = (ie-1)*ny + iy
"""
function make_joint_t_row(te::AbstractMatrix{<:Real}, ty::AbstractMatrix{<:Real})
    rowerr_e = maximum(abs.(sum(te, dims=2) .- 1.0))
    rowerr_y = maximum(abs.(sum(ty, dims=2) .- 1.0))
    @assert rowerr_e < 1e-10
    @assert rowerr_y < 1e-10

    # If z=(e,y) with e fastest, row-stoch joint is kron(te, ty) or kron(ty, te)
    # For iz = (ie-1)*ny + iy (y fastest), joint = kron(te, ty)
    return kron(Float64.(te), Float64.(ty))
end

function makegrids_disaster(ppp::Vector{Float64};
                            pi_annual::Float64,
                            kappa_e_log::Float64,
                            periods_per_year::Int=4)

    g0 = makegrids(ppp)

    @assert hasproperty(g0, :ex)
    @assert hasproperty(g0, :te)
    @assert hasproperty(g0, :ty)

    te0 = g0.te
    ty0 = g0.ty

    te_dis, meta = make_te_disaster_from_base_row(te0, g0.ex;
        pi_annual=pi_annual,
        periods_per_year=periods_per_year,
        kappa_e_log=kappa_e_log
    )

    t_dis = make_joint_t_row(te_dis, ty0)

    g_dis = (; g0..., te = te_dis, t = t_dis)
    return g_dis, meta
end
