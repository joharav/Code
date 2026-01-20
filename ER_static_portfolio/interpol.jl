# ==========================================================================
# 4D MODEL: Interpolation helpers (e,y discrete; w,d continuous)
# ==========================================================================

@inline function brack1d(grid::AbstractVector{<:Real}, x::Real)
    @inbounds begin
        n = length(grid)
        n == 0 && error("brack1d: empty grid")

        if x <= grid[1]
            return 1, 1, 0.0
        elseif x >= grid[n]
            return n, n, 0.0
        else
            j = searchsortedlast(grid, x)         # 1 <= j <= n-1
            gL = grid[j]
            gU = grid[j+1]
            w = (gU == gL) ? 0.0 : (x - gL) / (gU - gL)
            return j, j+1, w
        end
    end
end

# For discrete Markov states: map a value to the closest grid point index.
# If your simulation already has indices, you don't need this at all.
@inline function nearest_index(grid::AbstractVector{<:Real}, x::Real)
    @inbounds begin
        n = length(grid)
        n == 0 && error("nearest_index: empty grid")

        if x <= grid[1]
            return 1
        elseif x >= grid[n]
            return n
        else
            j = searchsortedlast(grid, x)     # 1..n-1
            # choose closest of j and j+1
            return (x - grid[j] <= grid[j+1] - x) ? j : (j+1)
        end
    end
end

# Main: given indices for (e,y), bilinear in (w,d)
@inline function interpol_ey(ie::Int, iy::Int, wold::Float64, dold::Float64,
                            g::NamedTuple, V::AbstractArray{<:Real,4})

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

# Optional convenience: map (e,y) values to indices then call interpol_ey.
# Use ONLY if you are not already tracking (ie,iy).
@inline function interpol(eold::Float64, yold::Float64, wold::Float64, dold::Float64,
                         g::NamedTuple, V::AbstractArray{<:Real,4})

    ie = nearest_index(g.ex, eold)
    iy = nearest_index(g.y,  yold)
    return interpol_ey(ie, iy, wold, dold, g, V)
end
