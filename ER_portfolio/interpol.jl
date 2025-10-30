# # interpol.jl — safe 5D (e,y,aa,a,d) multilinear interpolation for V (or any array)
# # Assumes g has fields: ex, y, aa, a, d (ascending, unique)
# # Returns a scalar interpolated value

# # --- helpers ---
# @inline clamp_to_grid(x::Float64, g::AbstractVector{<:Real}) = min(max(x, g[1]), g[end])

# @inline function _bracket(x::Float64, g::AbstractVector{<:Real})
#     @inbounds begin
#         n = length(g); @assert n ≥ 2 "Grid must have length ≥ 2"
#         if x ≤ g[1];        return 1, 2, 0.0; end
#         if x ≥ g[end];      return n-1, n, 1.0; end
#         j = searchsortedlast(g, x)
#         j = (j ≥ n ? n-1 : j)               # guard rare roundoff
#         w = (x - g[j]) / (g[j+1] - g[j])
#         return j, j+1, ifelse(w < 0, 0.0, ifelse(w > 1, 1.0, w))
#     end
# end

# # --- main ---
# function interpol(eold::Float64, yold::Float64, aaold::Float64, aold::Float64, dold::Float64,
#                   g::NamedTuple, v::Array{Float64,5})
#     @inbounds begin
#         # index ex,y on their own grids (if sim draws from those supports)
#         ie = max(1, min(searchsortedlast(g.ex, eold), length(g.ex)))
#         iy = max(1, min(searchsortedlast(g.y,  yold), length(g.y)))

#         # clamp continuous states to grids
#         aa = clamp_to_grid(aaold, g.aa)
#         a  = clamp_to_grid(aold,  g.a)
#         d  = clamp_to_grid(dold,  g.d)

#         # brackets + weights
#         iaL, iaU, wa    = _bracket(a,  g.a)
#         iAAL,iAAU,wAA   = _bracket(aa, g.aa)
#         idL, idU, wD    = _bracket(d,  g.d)

#         # 3D tri-linear on (aa, a, d), with fixed (ie,iy)
#         v000 = v[ie,iy,iAAL,iaL,idL]; v001 = v[ie,iy,iAAL,iaL,idU]
#         v010 = v[ie,iy,iAAL,iaU,idL]; v011 = v[ie,iy,iAAL,iaU,idU]
#         v100 = v[ie,iy,iAAU,iaL,idL]; v101 = v[ie,iy,iAAU,iaL,idU]
#         v110 = v[ie,iy,iAAU,iaU,idL]; v111 = v[ie,iy,iAAU,iaU,idU]

#         va0 = (1-wa)*v000 + wa*v010
#         va1 = (1-wa)*v001 + wa*v011
#         vb0 = (1-wa)*v100 + wa*v110
#         vb1 = (1-wa)*v101 + wa*v111

#         vd0 = (1-wAA)*va0 + wAA*vb0
#         vd1 = (1-wAA)*va1 + wAA*vb1

#         return (1-wD)*vd0 + wD*vd1
#     end
# end

function interpol(eold::Float64, yold::Float64, aaold::Float64, aold::Float64, dold::Float64,
    g::NamedTuple, v::Array{Float64,5})
    my_eps = 1.0e-6
    ide = findall(eold .== g.ex)
    idy = findall(yold .== g.y)

    # locate aa (local)
    iaalo = max(findall((aaold - my_eps) .<= g.aa)[1]-1, 1)
    if iaalo == sz.na
        iaahi = sz.na; aafrac = 1.0
    else
        iaahi = iaalo + 1
        aa_hi = g.aa[iaahi]; aa_lo = g.aa[iaalo]
        aafrac = abs(aaold - aa_lo) / abs(aa_hi - aa_lo)
    end

    # locate a (foreign)
    ialo = max(findall((aold - my_eps) .<= g.a)[1]-1, 1)
    if ialo == sz.na
        iahi = sz.na; afrac = 1.0
    else
        iahi = ialo + 1
        a_hi = g.a[iahi]; a_lo = g.a[ialo]
        afrac = abs(aold - a_lo) / abs(a_hi - a_lo)
    end

    # locate d
    iddlo = max(findall((dold - my_eps) .<= g.d)[1]-1, 1)
    if iddlo == sz.nd
        iddhi = sz.nd; dfrac = 1.0
    else
        iddhi = iddlo + 1
        d_hi = g.d[iddhi]; d_lo = g.d[iddlo]
        dfrac = abs(dold - d_lo) / abs(d_hi - d_lo)
    end

    # trilinear over (aa, a, d)
    v000 = v[ide, idy, iaalo, ialo, iddlo]
    v001 = v[ide, idy, iaalo, ialo, iddhi]
    v010 = v[ide, idy, iaalo, iahi, iddlo]
    v011 = v[ide, idy, iaalo, iahi, iddhi]
    v100 = v[ide, idy, iaahi, ialo, iddlo]
    v101 = v[ide, idy, iaahi, ialo, iddhi]
    v110 = v[ide, idy, iaahi, iahi, iddlo]
    v111 = v[ide, idy, iaahi, iahi, iddhi]

    va0 = (1-afrac)*v000 + afrac*v010
    va1 = (1-afrac)*v001 + afrac*v011
    vb0 = (1-afrac)*v100 + afrac*v110
    vb1 = (1-afrac)*v101 + afrac*v111

    v_d0 = (1-aafrac)*va0 + aafrac*vb0
    v_d1 = (1-aafrac)*va1 + aafrac*vb1

    vprime = (1-dfrac)*v_d0 + dfrac*v_d1
    return vprime
end
