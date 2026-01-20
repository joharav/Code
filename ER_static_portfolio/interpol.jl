# ==========================================================================
# 4D MODEL: Interpolation functions
# State: (e, y, w, d)
# ==========================================================================

# Bracket a value in a sorted grid
@inline function brack1d(grid::AbstractVector{<:Real}, x::Real)
    @inbounds begin
        n = length(grid)
        if n == 0
            error("brack1d: empty grid")
        end
        if x <= grid[1]
            return 1, 1, 0.0
        elseif x >= grid[n]
            return n, n, 0.0
        else
            j = searchsortedlast(grid, x)
            gL = grid[j]
            gU = grid[j+1]
            w = (gU == gL) ? 0.0 : (x - gL) / (gU - gL)
            return j, j+1, w
        end
    end
end

# Nearest index for discrete grids
@inline function near1d(grid::AbstractVector{<:Real}, x::Real)
    j = searchsortedlast(grid, x)
    j = ifelse(j < 1, 1, ifelse(j >= length(grid), length(grid), j))
    return j
end


# ==========================================================================
# Main interpolation function for 4D value/policy arrays
# V[ie, iy, iw, id]
# ==========================================================================
@inline function interpol(eold::Float64, yold::Float64, wold::Float64, dold::Float64,
                         g::NamedTuple, V::Array{Float64,4})
    
    # Discrete dimensions: nearest index
    ie = near1d(g.ex, eold)
    iy = near1d(g.y, yold)
    
    # Continuous dimensions: bracket and interpolate
    iwL, iwU, ww = brack1d(g.w, wold)
    idL, idU, wd = brack1d(g.d, dold)
    
    # Bilinear interpolation over (w, d)
    @inbounds begin
        v00 = V[ie, iy, iwL, idL]
        v10 = V[ie, iy, iwU, idL]
        v01 = V[ie, iy, iwL, idU]
        v11 = V[ie, iy, iwU, idU]
        
        v0 = (1.0 - ww) * v00 + ww * v10
        v1 = (1.0 - ww) * v01 + ww * v11
        
        return (1.0 - wd) * v0 + wd * v1
    end
end


# ==========================================================================
# Interpolation with separate e and y indices (for use in simulation)
# ==========================================================================
@inline function interpol_ey(ie::Int, iy::Int, wold::Float64, dold::Float64,
                            g::NamedTuple, V::Array{Float64,4})
    
    iwL, iwU, ww = brack1d(g.w, wold)
    idL, idU, wd = brack1d(g.d, dold)
    
    @inbounds begin
        v00 = V[ie, iy, iwL, idL]
        v10 = V[ie, iy, iwU, idL]
        v01 = V[ie, iy, iwL, idU]
        v11 = V[ie, iy, iwU, idU]
        
        v0 = (1.0 - ww) * v00 + ww * v10
        v1 = (1.0 - ww) * v01 + ww * v11
        
        return (1.0 - wd) * v0 + wd * v1
    end
end


# ==========================================================================
# Bracket helper (same as in maxbellman, but consolidated here)
# ==========================================================================
@inline function bracket_grid(x::Float64, grid::Vector{Float64})
    n = length(grid)
    @inbounds begin
        if x <= grid[1]
            return 1, 1, 0.0
        elseif x >= grid[n]
            return n, n, 0.0
        else
            j = searchsortedlast(grid, x)
            xL = grid[j]
            xU = grid[j+1]
            wt = (xU == xL) ? 0.0 : (x - xL) / (xU - xL)
            return j, j+1, wt
        end
    end
end
