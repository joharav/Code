# ==========================================================================
# 4D MODEL: Howard policy evaluation (ADJUSTMENT regime)
# Given a FIXED policy (w', d', s), update v under that policy.
#
# NOTE: queue must be E_{y'|y}[ V(e', y', w, d) ] already (your code does that).
# Here we only integrate over e' using trans_e.
# ==========================================================================

function howard(queue::Array{Float64,4}, util::Array{Float64,6},
    gidx::dtp.Ipol, grids::NamedTuple, pea::Vector{Float64})

beta    = pea[1]
rr      = (1 / beta) - 1
rr_star = pea[9]
kappa   = pea[11]

e_grid  = grids.ex
w_grid  = grids.w
wp_grid = grids.wp
s_grid  = grids.s
trans_e = grids.te

vnew = fill(-Inf, sz.ne, sz.ny, sz.nw, sz.nd)

Threads.@threads for id in 1:sz.nd
for iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
iwp = gidx.w[ie, iy, iw, id]
idp = gidx.d[ie, iy, iw, id]
is  = gidx.s[ie, iy, iw, id]

# Guard indices (prevents silent OOB corruption)
if !(1 <= iwp <= sz.npw && 1 <= idp <= sz.npd && 1 <= is <= sz.ns)
    continue
end

u_flow = util[ie, iy, iw, id, iwp, idp]
if u_flow <= -1e9
    vnew[ie, iy, iw, id] = -1e10
    continue
end

w_next = wp_grid[iwp]
s      = s_grid[is]
E_now  = e_grid[ie]

trans_cost = kappa * s * w_next

EV = 0.0
@inbounds for ie_next in 1:sz.ne
    E_next = e_grid[ie_next]

    w_realized = (1.0 - s) * w_next * (1.0 + rr) +
                 s * w_next * (1.0 + rr_star) * (E_next / E_now) -
                 trans_cost

    w_realized = clamp(w_realized, w_grid[1], w_grid[end])

    iw_L, iw_U, wt = bracket_grid(w_realized, w_grid)
    V_interp = (1.0 - wt) * queue[ie_next, iy, iw_L, idp] +
               wt        * queue[ie_next, iy, iw_U, idp]

    EV += trans_e[ie, ie_next] * V_interp
end

vnew[ie, iy, iw, id] = u_flow + beta * EV
end
end

return vnew, gidx
end