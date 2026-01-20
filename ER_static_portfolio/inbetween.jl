# Returns (lo, hi, w) such that:
#   x ≈ (1-w)*grid[lo] + w*grid[hi]
# with flat extrapolation at endpoints.
@inline function brack1d_idx(grid::AbstractVector{<:Real}, x::Real)
    @inbounds begin
        n = length(grid)
        n == 0 && error("brack1d_idx: empty grid")

        if x <= grid[1]
            return 1, 1, 0.0
        elseif x >= grid[n]
            return n, n, 0.0
        else
            j = searchsortedlast(grid, x)          # 1..n-1
            gL = grid[j]
            gU = grid[j+1]
            w  = (gU == gL) ? 0.0 : (x - gL) / (gU - gL)
            return j, j+1, w
        end
    end
end

function inbetween(grids::NamedTuple, islog::Bool, grid_type::Symbol)
    if grid_type == :w
        grid        = grids.w
        grid_policy = grids.wp
    elseif grid_type == :d
        grid        = grids.d
        grid_policy = grids.dp
    elseif grid_type == :a
        grid        = grids.a
        grid_policy = grids.ap
    elseif grid_type == :aa
        grid        = grids.aa
        grid_policy = grids.aap
    else
        error("Invalid grid type. Use :w, :d, :a, or :aa.")
    end

    if islog
        # log requires strictly positive support
        @assert minimum(grid) > 0 "inbetween: grid has nonpositive values, cannot log"
        @assert minimum(grid_policy) > 0 "inbetween: policy grid has nonpositive values, cannot log"
        grid        = log.(grid)
        grid_policy = log.(grid_policy)
    end

    np = length(grid_policy)
    lo = Vector{Int}(undef, np)
    hi = Vector{Int}(undef, np)
    w  = Vector{Float64}(undef, np)

    @inbounds for i in 1:np
        l, h, wi = brack1d_idx(grid, grid_policy[i])
        lo[i] = l
        hi[i] = h
        w[i]  = wi
    end

    return (lo = lo, hi = hi, w = w)
end
