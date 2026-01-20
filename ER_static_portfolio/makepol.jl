# ==========================================================================
# 4D MODEL: Policy extraction functions
# ==========================================================================

# Convert policy indices to levels
function makepol(gidx::Array{Int,4}, pgrid::Vector{Float64})
    pol = zeros(size(gidx))
    @inbounds for i in eachindex(gidx)
        pol[i] = pgrid[gidx[i]]
    end
    return pol
end


# ==========================================================================
# Compute consumption from (w, d, s) policies
# adjust = 1 for adjustment regime, 0 for non-adjust
# ==========================================================================
function makepol_c(pol_w::Array{Float64,4}, pol_d::Array{Float64,4}, 
                  pol_s::Array{Float64,4}, grids::NamedTuple, 
                  adjust::Int, pea::Vector{Float64})
    
    delta = pea[2]
    f = pea[7]
    wage = pea[8]
    pd = pea[10]
    kappa = pea[11]
    tau = pea[12]
    h = pea[13]
    chi = pea[16]
    ft = pea[17]
    
    w_grid = grids.w
    d_grid = grids.d
    e_grid = grids.ex
    y_grid = grids.y
    
    c = zeros(sz.ne, sz.ny, sz.nw, sz.nd)
    
    Threads.@threads for id in 1:sz.nd
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    E = e_grid[ie]
                    Y = y_grid[iy]
                    w_now = w_grid[iw]
                    d_now = d_grid[id]
                    
                    w_next = pol_w[ie, iy, iw, id]
                    d_next = pol_d[ie, iy, iw, id]
                    s_next = pol_s[ie, iy, iw, id]
                    
                    # Income
                    income = Y * wage * h * (1 - tau) + w_now
                    
                    if adjust == 1
                        # Adjustment regime: sell old durables, buy new
                        sale_value = E * pd * (1 - f) * (1 - delta) * d_now
                        durable_purchase = E * pd * d_next
                        time_cost = wage * h * ft * Y
                        
                        c[ie, iy, iw, id] = income + sale_value - durable_purchase - 
                                           w_next - time_cost
                    else
                        # Non-adjust: only maintenance
                        maintenance = E * pd * delta * chi * d_now
                        c[ie, iy, iw, id] = income - w_next - maintenance
                    end
                end
            end
        end
    end
    
    return c
end


# ==========================================================================
# Derive peso and dollar asset holdings from total wealth and dollar share
# ==========================================================================
function derive_asset_holdings(pol_w::Array{Float64,4}, pol_s::Array{Float64,4}, 
                              grids::NamedTuple)
    
    e_grid = grids.ex
    
    aa = zeros(sz.ne, sz.ny, sz.nw, sz.nd)  # peso assets
    a = zeros(sz.ne, sz.ny, sz.nw, sz.nd)   # dollar assets
    
    for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
        w_pr = pol_w[ie, iy, iw, id]
        s_pr = pol_s[ie, iy, iw, id]
        E = e_grid[ie]
        
        # w' = aa' + E * a'
        # s = E * a' / w'  =>  a' = s * w' / E
        # aa' = (1 - s) * w'
        
        aa[ie, iy, iw, id] = (1.0 - s_pr) * w_pr
        a[ie, iy, iw, id] = (s_pr * w_pr) / max(E, 1e-10)
    end
    
    return aa, a
end
