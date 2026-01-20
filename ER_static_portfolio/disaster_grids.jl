module DisasterGrids

using Distributions, LinearAlgebra

# Tauchen in logs with optional jump shift
function tauchen_log(ne::Int, ρ::Float64, σ::Float64;
                     nstd::Float64=3.0,
                     μ_ss::Float64=0.0,
                     x_grid::Union{Nothing,Vector{Float64}}=nothing,
                     jump::Float64=0.0)

    σx = σ / sqrt(1 - ρ^2)

    if x_grid === nothing
        x_min = μ_ss - nstd * σx
        x_max = μ_ss + nstd * σx
        x_grid = range(x_min, x_max; length=ne) |> collect
    else
        @assert length(x_grid) == ne "x_grid length != ne in tauchen_log"
    end

    dx = x_grid[2] - x_grid[1]
    edges = similar(x_grid, ne+1)
    edges[1]   = x_grid[1] - dx/2
    edges[end] = x_grid[end] + dx/2
    @inbounds for i in 2:ne
        edges[i] = (x_grid[i-1] + x_grid[i]) / 2
    end

    P   = zeros(Float64, ne, ne)
    N01 = Normal(0.0, 1.0)

    @inbounds for i in 1:ne
        μ = ρ * x_grid[i] + jump
        for j in 1:ne
            z_hi = (edges[j+1] - μ) / σ
            z_lo = (edges[j]   - μ) / σ
            P[i,j] = cdf(N01, z_hi) - cdf(N01, z_lo)
        end
        rs = sum(P[i,:])
        if rs <= 0
            P[i,:] .= 1.0/ne
        else
            P[i,:] ./= rs
        end
    end

    return x_grid, P
end

# annual → quarterly disaster prob
pi_quarterly(pi_y::Float64) = 1.0 - (1.0 - pi_y)^(0.25)

"""
    build_ex_grid_disaster(ρ_e, σ_e; ne, nstd, pi_annual, kappa_e, μ_ss_loge=0.0)

log e' = ρ_e log e + σ_e ε + J κ_e, J ~ Bernoulli(π_q),
implemented as P = (1-π_q) P_base + π_q P_jump.
"""
function build_ex_grid_disaster(ρ_e::Float64,
                                σ_e::Float64;
                                ne::Int,
                                nstd::Float64,
                                pi_annual::Float64,
                                kappa_e::Float64,
                                μ_ss_loge::Float64=0.0)

    πq = pi_quarterly(pi_annual)

    # baseline
    x_grid, P_base = tauchen_log(ne, ρ_e, σ_e;
                                 nstd=nstd,
                                 μ_ss=μ_ss_loge,
                                 x_grid=nothing,
                                 jump=0.0)

    # jump case on same grid
    _, P_jump = tauchen_log(ne, ρ_e, σ_e;
                            nstd=nstd,
                            μ_ss=μ_ss_loge,
                            x_grid=x_grid,
                            jump=kappa_e)

    P_ex = (1.0 - πq) .* P_base .+ πq .* P_jump

    @inbounds for i in 1:ne
        s = sum(P_ex[i,:])
        if s <= 0
            P_ex[i,:] .= 1.0/ne
        else
            P_ex[i,:] ./= s
        end
    end

    ex_grid = exp.(x_grid)
    return ex_grid, P_ex
end

end # module

using .DisasterGrids

"""
    makegrids_disaster(ppp; pi_annual, kappa_e)

Same as `makegrids(ppp)` but with an exchange-rate rare-disaster process:
log e' = ρ_e log e + σ_e ε + J κ_e, J ~ Bernoulli(π_e).
"""
function makegrids_disaster(ppp::Vector{Float64};
                            pi_annual::Float64,
                            kappa_e::Float64)

    # 1. Baseline grids (assets, durables, etc.)
    g0 = makegrids(ppp)

    # 2. Pull parameters needed for e,y processes
    rho_e   = ppp[3]
    sigma_e = ppp[4]
    rho_y   = ppp[14]
    sigma_y = ppp[15]

    # 3. Disaster exchange-rate process (on same dimension ne)
    @assert sz.ne > 2 "Disaster grid only coded for sz.ne > 2"
    ex_grid, P_ex = DisasterGrids.build_ex_grid_disaster(
        rho_e, sigma_e;
        ne        = sz.ne,
        nstd      = sz.nstd_e,
        pi_annual = pi_annual,
        kappa_e   = kappa_e,
        μ_ss_loge = 3.0,
    )

    # 4. Income process: replicate exactly your makegrids logic
    if sz.ny == 2
        trans_y = [0.9 0.1; 0.8 0.2]
        yg      = zeros(2)
        yg[1]   = 0.90
        yg[2]   = 1.10
    else
        numy     = sz.ny
        numstd_y = sz.nstd_y
        mew      = 0.0
        yg, trans_y = tauchen(mew, sigma_y, rho_y, numy, numstd_y)
        yg = exp.(yg)
    end

    # 5. Joint (e,y) transition under disaster
    trans_disaster = kron(P_ex, trans_y)

    # 6. Return same structure as g0, but with new ex, y, t
    return (; g0..., ex = ex_grid, y = yg, t = trans_disaster)
end
