# ---------- tiny helpers ----------
@inline clamp_to_grid(x, g::AbstractVector{<:Real}) =
    min(max(x, first(g)), last(g))

# Sample an index from a CDF row (monotone, ends at 1). Returns 1..length(cdfrow)
@inline function draw_next(cdfrow::AbstractVector{<:Real}, u::Real)
    @inbounds begin
        n  = length(cdfrow)
        @assert n > 0 "draw_next: empty cdf row"
        uu = min(max(u, 0.0), 1.0 - eps())   # keep in [0,1)
        j  = searchsortedfirst(cdfrow, uu)   # first index with cdf ≥ uu
        return j < 1 ? 1 : (j > n ? n : j)
    end
end

function simmodel(answ::NamedTuple)
    # Unpack once
    grids   = answ.g
    tmat    = grids.t
    exg     = grids.ex
    yg      = grids.y

    # Precompute joint (e,y) CDF rows for transitions
    phatcdf = cumsum(tmat, dims=2)           # size: (ne*ny, ne*ny)
    @inbounds phatcdf[:, end] .= 1.0

    # Crude ergodic-ish initializer; force well-formed CDF
    cdf_wgt = Matrix(tmat')^100
    cdf_wgt = cumsum(cdf_wgt[:, Int(floor(sz.ne * 0.5)) + 1])
    @inbounds cdf_wgt[end] = 1.0

    # Storage
    allv        = zeros(sz.nYears, sz.nFirms)
    alla        = zeros(sz.nYears, sz.nFirms)     # foreign asset
    allaa       = zeros(sz.nYears, sz.nFirms)     # local asset
    alle        = zeros(sz.nYears, sz.nFirms)
    ally        = zeros(sz.nYears, sz.nFirms)
    alld        = zeros(sz.nYears, sz.nFirms)
    allc        = zeros(sz.nYears, sz.nFirms)
    adjust_indicator = zeros(sz.nYears, sz.nFirms)

    # Initial joint-state indices (1..ne*ny), one per firm
    ls = zeros(Int, sz.nYears + 1, sz.nFirms)
    @inbounds for ifi in 1:sz.nFirms
        ls[1, ifi] = draw_next(cdf_wgt, globals.draws[1, ifi])
    end

    # Initial idiosyncratic indices on grids
    astart  = globals.draws[1, :]
    aastart = globals.draws[3, :]
    dstart  = globals.draws[2, :]


    Threads.@threads for ifi in 1:sz.nFirms
        @inbounds begin
            picka  = min(Int(floor(sz.na * astart[ifi]))  + 1, sz.na)
            pickaa = min(Int(floor(sz.na * aastart[ifi])) + 1, sz.na)
            pickd  = min(Int(floor(sz.nd * dstart[ifi]))  + 1, sz.nd)

            pickey = ls[1, ifi]
            picke  = div(pickey - 1, sz.ny) + 1
            picky  = mod(pickey - 1, sz.ny) + 1

            # Start from the merged-policy lattice point (only to get initial levels)
            a_old   = answ.pol.a[ picke, picky, pickaa, picka, pickd ]
            aa_old  = answ.pol.aa[picke, picky, pickaa, picka, pickd ]
            d_old   = answ.pol.d[ picke, picky, pickaa, picka, pickd ]
            v_old   = answ.v[     picke, picky, pickaa, picka, pickd ]
            c_old   = answ.pol.c[ picke, picky, pickaa, picka, pickd ]

            for iti in 1:sz.nYears
                # Aggregate states (discrete)
                e = exg[picke]
                y = yg[picky]

                # Clamp continuous states to grids before interpolation
                e_in  = e
                y_in  = clamp_to_grid(y,   grids.y)
                aa_in = clamp_to_grid(aa_old, grids.aa)
                a_in  = clamp_to_grid(a_old,  grids.a)
                d_in  = clamp_to_grid(d_old,  grids.d)

                # Compare values across regimes at the current state
                vA = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.adjust_result.v)
                vN = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.noadjust_result.v)
                do_adjust = (vA - vN) > sz.toler

                # Interpolate policies from the winner
                if do_adjust
                    a_pr  = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.adjust_result.pol.a)
                    aa_pr = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.adjust_result.pol.aa)
                    d_pr  = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.adjust_result.pol.d)
                    c_pr  = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.adjust_result.pol.c)
                    v_pr  = vA
                else
                    a_pr  = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.noadjust_result.pol.a)
                    aa_pr = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.noadjust_result.pol.aa)
                    d_pr  = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.noadjust_result.pol.d)
                    c_pr  = interpol(e_in, y_in, aa_in, a_in, d_in, grids, answ.noadjust_result.pol.c)
                    v_pr  = vN
                end

                # Record adjustment only if chosen AND the durable truly moves
                adjust_indicator[iti, ifi] = (do_adjust && abs(d_pr - d_in) > sz.toler) ? 1.0 : 0.0

                # Save paths
                allv[iti, ifi]  = v_pr
                alla[iti, ifi]  = a_pr
                allaa[iti, ifi] = aa_pr
                alle[iti, ifi]  = e
                ally[iti, ifi]  = y
                alld[iti, ifi]  = d_pr
                allc[iti, ifi]  = c_pr

                # Next (e,y)
                row = ls[iti, ifi]
                u   = globals.draws[iti + 1, ifi]
                nxt = draw_next(view(phatcdf, row, :), u)
                ls[iti + 1, ifi] = nxt
                picke = div(nxt - 1, sz.ny) + 1
                picky = mod(nxt - 1, sz.ny) + 1

                # Update continuous states
                v_old  = v_pr
                a_old  = a_pr
                aa_old = aa_pr
                d_old  = d_pr
                c_old  = c_pr
            end
        end
    end

    return (
        v  = allv,
        d  = alld,
        a  = alla,
        aa = allaa,
        ex = alle,
        y  = ally,
        c  = allc,
        adjust_indicator = adjust_indicator,
    )
end
