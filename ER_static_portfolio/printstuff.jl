using Printf

function printstuff(answ::NamedTuple)
    # --- unpack ---
    g   = answ.g
    outdir = "Output/Policy"
    isdir(outdir) || mkpath(outdir)

    # Value functions by regime (recommended; your answ.v is ambiguous now)
    vA = answ.adjust_result.v          # (ne,ny,nw,nd)
    vN = answ.noadjust_result.v        # (ne,ny,nw,nd)

    # Policies (4D)
    polA = answ.adjust_result.pol
    polN = answ.noadjust_result.pol

    # ---------- value functions ----------
    open(joinpath(outdir, "v_adjust.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", vA[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "v_noadjust.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", vN[ie, iy, iw, id])
        end
    end

    # ---------- policies (adjust) ----------
    open(joinpath(outdir, "pol_adjust_w.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polA.w[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "pol_adjust_d.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polA.d[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "pol_adjust_s.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polA.s[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "pol_adjust_c.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polA.c[ie, iy, iw, id])
        end
    end

    # ---------- policies (no-adjust) ----------
    open(joinpath(outdir, "pol_noadjust_w.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polN.w[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "pol_noadjust_d.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polN.d[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "pol_noadjust_s.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polN.s[ie, iy, iw, id])
        end
    end
    open(joinpath(outdir, "pol_noadjust_c.txt"), "w") do io
        for id in 1:sz.nd, iw in 1:sz.nw, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", polN.c[ie, iy, iw, id])
        end
    end

    # ---------- grids & transition ----------
    eg = g.ex
    yg = g.y
    wg = g.w
    dg = g.d
    sg = hasproperty(g, :s)  ? g.s  : Float64[]
    wpg = hasproperty(g, :wp) ? g.wp : Float64[]
    trans = g.t

    open(joinpath(outdir, "statespace.txt"), "w") do io
        @printf(io, "wealth state grid (w)\n")
        for x in wg; @printf(io, "%16.8f \n", x); end

        if !isempty(wpg)
            @printf(io, "\nwealth policy grid (wp)\n")
            for x in wpg; @printf(io, "%16.8f \n", x); end
        end

        @printf(io, "\ndurable state grid (d)\n")
        for x in dg; @printf(io, "%16.8f \n", x); end

        if !isempty(sg)
            @printf(io, "\nportfolio share grid (s)\n")
            for x in sg; @printf(io, "%16.8f \n", x); end
        end

        @printf(io, "\nexchange rate grid (e)\n")
        for x in eg; @printf(io, "%16.8f \n", x); end

        @printf(io, "\nincome grid (y)\n")
        for x in yg; @printf(io, "%16.8f \n", x); end

        @printf(io, "\ntransition matrix (z' × z), z=(e,y)\n")
        for jj in 1:(sz.ne*sz.ny)
            for ii in 1:(sz.ne*sz.ny)
                @printf(io, "%16.8f", trans[ii, jj])
            end
            @printf(io, "\n")
        end
    end

    return 1
end
