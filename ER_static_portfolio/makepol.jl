# ==========================================================================
# 4D MODEL: Policy extraction functions (CORRECT + CONSISTENT)
# State: (e, y, w, d)
# Policies stored as integer indices in gidx (Ipol), converted to levels here.
# ==========================================================================
using Base.Threads

# --------------------------------------------------------------------------
# Convert integer index policy to levels on its grid
# --------------------------------------------------------------------------
function makepol(gidx::AbstractArray{<:Integer,4}, pgrid::AbstractVector{<:Real})
    pol = Array{Float64}(undef, size(gidx))
    @inbounds for i in eachindex(gidx)
        pol[i] = Float64(pgrid[gidx[i]])
    end
    return pol
end



# --------------------------------------------------------------------------
# Derive peso (aa) and dollar (a) asset holdings from total wealth and dollar share
#
# Convention:
#   w' = aa' + E * a'
#   s  = (E * a') / w'
# so:
#   a'  = s*w'/E
#   aa' = (1-s)*w'
# --------------------------------------------------------------------------
function derive_asset_holdings(pol_w::Array{Float64,4}, pol_s::Array{Float64,4},
                               grids::NamedTuple)

    e_grid = grids.ex

    aa = zeros(Float64, sz.ne, sz.ny, sz.nw, sz.nd)
    a  = zeros(Float64, sz.ne, sz.ny, sz.nw, sz.nd)

    Threads.@threads for id in 1:sz.nd
        for iw in 1:sz.nw
            for iy in 1:sz.ny
                for ie in 1:sz.ne
                    w_pr = pol_w[ie, iy, iw, id]
                    s_pr = pol_s[ie, iy, iw, id]
                    E    = e_grid[ie]

                    aa[ie, iy, iw, id] = (1.0 - s_pr) * w_pr
                    a[ie, iy, iw, id]  = (s_pr * w_pr) / max(E, 1e-10)
                end
            end
        end
    end

    return aa, a
end
