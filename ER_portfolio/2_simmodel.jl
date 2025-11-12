# ---------- tiny helpers ----------
# clamp x to a grid's min/max
@inline clamp_to_grid(x, g::AbstractVector{<:Real}) =
    min(max(x, first(g)), last(g))

# sample an index from a CDF row (monotone, ends at 1). Always returns 1..length(cdfrow)
@inline function draw_next(cdfrow::AbstractVector{<:Real}, u::Real)
    @inbounds begin
        n  = length(cdfrow)
        @assert n > 0 "draw_next: empty cdf row"
        uu = min(max(u, 0.0), 1.0 - eps())   # keep in [0, 1)
        j  = searchsortedfirst(cdfrow, uu)   # first index with cdf ≥ uu
        return j < 1 ? 1 : (j > n ? n : j)
    end
end


function simmodel(answ::NamedTuple)
    # Unpack
    v           = answ.v
    pol         = answ.pol                 # expects fields: a, aa, d, c
    grids       = answ.g
    tmat        = grids.t
    d_adjust    = answ.adjust_result.pol.d
    adjflag_grid = answ.adjust_flag


    exg = grids.ex
    yg  = grids.y

    # Precompute joint (e,y) CDF rows for transitions
    phatcdf = cumsum(tmat, dims=2)                 # size: (ne*ny, ne*ny)
    @inbounds phatcdf[:, end] .= 1.0               # force exact 1.0 at the end

    # Optional: a crude stationary-ish init weight; keep but make safe
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0

    # Storage
    allv        = zeros(sz.nYears, sz.nFirms)
    alla        = zeros(sz.nYears, sz.nFirms)     # foreign
    allaa       = zeros(sz.nYears, sz.nFirms)     # local
    alle        = zeros(sz.nYears, sz.nFirms)
    ally        = zeros(sz.nYears, sz.nFirms)
    alld        = zeros(sz.nYears, sz.nFirms)
    alld_adjust = zeros(sz.nYears, sz.nFirms)
    allc        = zeros(sz.nYears, sz.nFirms)
    adjust_indicator        = zeros(sz.nYears, sz.nFirms)


    # Initial joint-state indices (1..ne*ny), one per firm
    ls = zeros(Int, sz.nYears + 1, sz.nFirms)
    @inbounds for ifi in 1:sz.nFirms
        u0 = globals.draws[1, ifi]
        ls[1, ifi] = draw_next(cdf_wgt, u0)
    end

    # initial idiosyncratic indices on asset/durable grids
    astart  = globals.draws[1, :]
    aastart = globals.draws[3, :]              # NEW: local asset init
    dstart  = globals.draws[2, :]

    Threads.@threads for ifi in 1:sz.nFirms
        @inbounds begin
            picka  = min(Int(floor(sz.na * astart[ifi]))  + 1, sz.na)
            pickaa = min(Int(floor(sz.na * aastart[ifi])) + 1, sz.na)
            pickd  = min(Int(floor(sz.nd * dstart[ifi]))  + 1, sz.nd)

            # unpack joint (e,y) from flat index
            pickey = ls[1, ifi]
            picke  = div(pickey - 1, sz.ny) + 1
            picky  = mod(pickey - 1, sz.ny) + 1

            # Start from policy at that lattice point
            aold  = pol.a[ picke, picky, pickaa, picka, pickd ]
            aaold = pol.aa[picke, picky, pickaa, picka, pickd ]
            dold  = pol.d[ picke, picky, pickaa, picka, pickd ]
            d_adjustold = d_adjust[picke, picky, pickaa, picka, pickd ]
            vold  = v[ picke, picky, pickaa, picka, pickd ]
            cold  = pol.c[picke, picky, pickaa, picka, pickd ]
            adj_old = answ.adjust_flag[picke, picky, pickaa, picka, pickd ]

            for iti in 1:sz.nYears
                eold = exg[picke]
                yold = yg[picky]

                # Clamp continuous states before interpolating
                e_in  = eold                                # if discrete, OK
                y_in  = clamp_to_grid(yold,  grids.y)
                aa_in = clamp_to_grid(aaold, grids.aa)
                a_in  = clamp_to_grid(aold,  grids.a)
                d_in  = clamp_to_grid(dold,  grids.d)

                # Interpolate continuation values and next policies (scalars)
                vprime        = interpol(e_in, y_in, aa_in, a_in, d_in, grids, v)
                aprime        = interpol(e_in, y_in, aa_in, a_in, d_in, grids, pol.a)
                aaprime       = interpol(e_in, y_in, aa_in, a_in, d_in, grids, pol.aa)
                dprime        = interpol(e_in, y_in, aa_in, a_in, d_in, grids, pol.d)
                d_adjustprime = interpol(e_in, y_in, aa_in, a_in, d_in, grids, d_adjust)
                cprime        = interpol(e_in, y_in, aa_in, a_in, d_in, grids, pol.c)
                adjflag_prime = interpol(e_in, y_in, aa_in, a_in, d_in, grids, adjflag_grid)



                # Draw next joint (e,y) using the row CDF
                row = ls[iti, ifi]
                u   = globals.draws[iti + 1, ifi]
                nxt = draw_next(view(phatcdf, row, :), u)
                ls[iti + 1, ifi] = nxt
                picke = div(nxt - 1, sz.ny) + 1
                picky = mod(nxt - 1, sz.ny) + 1

                eprime = exg[picke]
                yprime = yg[picky]

                # Save
                allv[iti, ifi]        = vprime
                alla[iti, ifi]        = aprime
                allaa[iti, ifi]       = aaprime
                alle[iti, ifi]        = eprime
                ally[iti, ifi]        = yprime
                alld[iti, ifi]        = dprime
                alld_adjust[iti, ifi] = d_adjustprime
                allc[iti, ifi]        = cprime
                adjust_indicator[iti, ifi] = adjflag_prime ≥ 0.5 ? 1.0 : 0.0


                # Update state
                vold  = vprime
                aold  = aprime
                aaold = aaprime
                dold  = dprime
                d_adjustold = d_adjustprime
                cold  = cprime
                adj_old = adjflag_prime

            end
        end
    end


    return (
        v = allv,
        d = alld,
        a = alla,
        aa = allaa,
        ex = alle,
        y = ally,
        c = allc,
        adjust_indicator = adjust_indicator
    )
end
