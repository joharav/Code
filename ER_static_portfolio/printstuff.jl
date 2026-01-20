function printstuff(answ::NamedTuple)
    vnew  = answ.v
    pol   = answ.pol
    apol  = pol.a
    aapol = pol.aa
    dpol  = pol.d
    cpol  = pol.c
    g     = answ.g

    # ensure dir
    outdir = "Output/Policy"
    isdir(outdir) || mkpath(outdir)

    # ---------- value function ----------
    open(joinpath(outdir, "v.txt"), "w") do io
        for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", vnew[ie, iy, iaa, ia, id])
        end
    end

    # ---------- policies ----------
    open(joinpath(outdir, "a.txt"), "w") do io
        for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", apol[ie, iy, iaa, ia, id])
        end
    end
    open(joinpath(outdir, "aa.txt"), "w") do io
        for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", aapol[ie, iy, iaa, ia, id])
        end
    end
    open(joinpath(outdir, "d.txt"), "w") do io
        for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", dpol[ie, iy, iaa, ia, id])
        end
    end
    open(joinpath(outdir, "c.txt"), "w") do io
        for id in 1:sz.nd, ia in 1:sz.na, iaa in 1:sz.na, iy in 1:sz.ny, ie in 1:sz.ne
            @printf(io, "%16.8f\n", cpol[ie, iy, iaa, ia, id])
        end
    end

    # ---------- grids & transition ----------
    eg, yg   = g.ex, g.y
    ag, apg  = g.a, g.ap
    aag, aap = g.aa, g.aap
    dg, dpg  = g.d, g.dp
    trans    = g.t

    open(joinpath(outdir, "statespace.txt"), "w") do io
        @printf(io, "foreign asset state grid (a)\n")
        for x in ag;  @printf(io, "%16.8f \n", x); end
        @printf(io, "\nforeign asset policy grid (ap)\n")
        for x in apg; @printf(io, "%16.8f \n", x); end

        @printf(io, "\nlocal asset state grid (aa)\n")
        for x in aag; @printf(io, "%16.8f \n", x); end
        @printf(io, "\nlocal asset policy grid (aap)\n")
        for x in aap; @printf(io, "%16.8f \n", x); end

        @printf(io, "\ndurable state grid\n")
        for x in dg;  @printf(io, "%16.8f \n", x); end
        @printf(io, "\ndurable policy grid\n")
        for x in dpg; @printf(io, "%16.8f \n", x); end

        @printf(io, "\nexchange rate grid\n")
        for x in eg;  @printf(io, "%16.8f \n", x); end
        @printf(io, "\nidiosyncratic income grid\n")
        for x in yg;  @printf(io, "%16.8f \n", x); end

        @printf(io, "\ntransition matrix (z' × z)\n")
        for jj in 1:(sz.ne*sz.ny)
            for ii in 1:(sz.ne*sz.ny)
                @printf(io, "%16.8f", trans[ii, jj])
            end
            @printf(io, "\n")
        end
    end

    return 1
end
