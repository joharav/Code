# --------------------------------------------------------------------------
# Consumption implied by policies
#
# IMPORTANT accounting convention (matches your utility + simmodel):
# - State variable w is total liquid wealth *available this period* (already
#   incorporating last period returns and FX valuation).
# - Choice w' is next-period total liquid wealth saved at end of period.
# - Portfolio share s affects NEXT period w via returns/FX, NOT today's c.
#
# adjust = 1 (adjust regime): sell old durable, buy new, pay time cost
# adjust = 0 (no-adjust): keep durable stock with maintenance cost
#
# NOTE: kappa not in c, because you already charge transaction costs in wealth
#       transition in Bellman/simulation (w_{t+1} law of motion).
# --------------------------------------------------------------------------

function makepol_c(pol_w::Array{Float64,4}, pol_d::Array{Float64,4},
                   grids::NamedTuple, adjust::Int, pea::Vector{Float64})

    delta = pea[2]
    f     = pea[7]
    wage  = pea[8]
    pd    = pea[10]
    tau   = pea[12]
    h     = pea[13]
    chi   = pea[16]
    ft    = pea[17]

    w_grid = grids.w
    d_grid = grids.d
    e_grid = grids.ex
    y_grid = grids.y

    c = Array{Float64}(undef, sz.ne, sz.ny, sz.nw, sz.nd)

    Threads.@threads for id in 1:sz.nd
        d_now = d_grid[id]
        for iw in 1:sz.nw
            w_now = w_grid[iw]
            for iy in 1:sz.ny
                Y = y_grid[iy]
                labor_inc = Y * wage * h * (1 - tau)
                for ie in 1:sz.ne
                    E = e_grid[ie]

                    w_next = pol_w[ie, iy, iw, id]
                    d_next = pol_d[ie, iy, iw, id]

                    if adjust == 1
                        sale_value       = E * pd * (1 - f) * (1 - delta) * d_now
                        durable_purchase = E * pd * d_next
                        time_cost        = wage * h * ft * Y
                        c[ie, iy, iw, id] = labor_inc + w_now + sale_value - durable_purchase - w_next - time_cost
                    else
                        maintenance      = E * pd * delta * chi * d_now
                        c[ie, iy, iw, id] = labor_inc + w_now - w_next - maintenance
                    end
                end
            end
        end
    end

    return c
end