# === interpol_helpers.jl ===
# Deterministic bracketing with clamping. Never returns empty.
@inline function brack1d(grid::AbstractVector{<:Real}, x::Real)
    @inbounds begin
        n = length(grid)
        @assert n > 0 "brack1d: empty grid"
        if x ≤ grid[1]             # left clamp
            return 1, 1, 0.0
        elseif x ≥ grid[n]         # right clamp
            return n, n, 0.0
        else
            j = searchsortedlast(grid, x)  # 1 ≤ j ≤ n-1 even with ties
            gL = grid[j]; gU = grid[j+1]
            w  = (gU == gL) ? 0.0 : (x - gL) / (gU - gL)
            return j, j+1, w
        end
    end
end

# For discrete grids where you just want the nearest index (e.g., shocks)
@inline function near1d(grid::AbstractVector{<:Real}, x::Real)
    j = searchsortedlast(grid, x)
    j = ifelse(j < 1, 1, ifelse(j ≥ length(grid), length(grid), j))
    # pick j or j+1 — here we stick with j (left) to keep it deterministic
    return j
end
# === interpol.jl ===
# Old (problematic): using findall / count / idx[1]
# idx = findall(g.aa .≤ aaold .≤ g.aa[end])  # etc...
# i = idx[1]   # <-- blows up when idx == []

# New (safe): use brack1d for continuous, near1d for discrete
@inline function interpol(eold::Float64, yold::Float64, aaold::Float64,
    aold::Float64,  dold::Float64,
    g::NamedTuple,   V::Array{Float64,5})

# Discrete (likely Markov) dimensions: pick nearest *valid* index
ie = near1d(g.ex, eold)
iy = near1d(g.y,  yold)

# Continuous dims: bracket + weight (always valid)
iaaL, iaaU, waa = brack1d(g.aa, aaold)
iaL,  iaU,  wa  = brack1d(g.a,  aold)
idL,  idU,  wd  = brack1d(g.d,  dold)

# Trilinear blend over (aa, a, d), holding (e,y) fixed
@inbounds begin
v000 = V[ie,iy,iaaL,iaL, idL]
v100 = V[ie,iy,iaaU,iaL, idL]
v010 = V[ie,iy,iaaL,iaU, idL]
v110 = V[ie,iy,iaaU,iaU, idL]

v001 = V[ie,iy,iaaL,iaL, idU]
v101 = V[ie,iy,iaaU,iaL, idU]
v011 = V[ie,iy,iaaL,iaU, idU]
v111 = V[ie,iy,iaaU,iaU, idU]

v0 = (1-waa)*( (1-wa)*v000 + wa*v010 ) + waa*( (1-wa)*v100 + wa*v110 )
v1 = (1-waa)*( (1-wa)*v001 + wa*v011 ) + waa*( (1-wa)*v101 + wa*v111 )
return (1-wd)*v0 + wd*v1
end
end
