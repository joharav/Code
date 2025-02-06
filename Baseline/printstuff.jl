function printstuff(answ::NamedTuple)

    vnew = answ.v

    pol = answ.pol
    apol = pol.a
    dpol = pol.d
    grids = answ.g

    filename = "Baseline/Output/v.txt"
    io = open(filename, "w")

    for id in 1:sz.nd
        for ia in 1:sz.na
            for ip in 1:sz.np
                @printf(io, "%16.8f", vnew[ip, ia, id])
                if ip < sz.np
                    @printf(io, " ")  # Space between p values
                end
            end
            @printf(io, "\n")  # New line for each d within an a block
        end
        @printf(io, "\n")  # Extra newline to separate blocks of a for readability
    end

    close(io)


    filname = "Baseline/Output/a.txt"
    io = open(filname, "w")

    for id in 1:sz.nd
        for ia in 1:sz.na
            for ip in 1:sz.np
                @printf(io, "%16.8f", apol[ip, ia, id])
                #if iz < sz.nz
                #    @printf(io, " ")  # Space between z values
                #end
            end
            @printf(io, "\n")  # New line for each d within an a block
        end
        @printf(io, "\n")  # Extra newline to separate blocks of a for readability
    end
    close(io)


    filname = "Baseline/Output/d.txt"
    io = open(filname, "w")
    for id in 1:sz.nd
        for ia in 1:sz.na
            for ip in 1:sz.np
                @printf(io, "%16.8f", dpol[ip, ia, id])
                #if iz < sz.nz
                #    @printf(io, " ")  # Space between p values
                #end
            end
            @printf(io, "\n")  # New line for each d within an a block
        end
        @printf(io, "\n")  # Extra newline to separate blocks of a for readability
    end
    close(io)

    pg = grids.p;
    ag = grids.a;
    apg = grids.ap;
    dg = grids.d;
    dpg = grids.dp;
    trans = grids.t; 
    filname = "Baseline/Output/statespace.txt"
    io = open(filname,"w")
    @printf(io," asset state grid\n")
    for jj in 1:sz.na
        @printf(io,"%16.8f \n",ag[jj])  
    end
    @printf(io," \n")
    @printf(io," asset policy grid\n")
    for jj in 1:sz.npa
        @printf(io,"%16.8f \n",apg[jj])  
    end
    @printf(io," \n")
    @printf(io," durable state grid\n")
    for jj in 1:sz.na
        @printf(io,"%16.8f \n",dg[jj])  
    end
    @printf(io," \n")
    @printf(io," durable policy grid\n")
    for jj in 1:sz.npa
        @printf(io,"%16.8f \n",dpg[jj])  
    end
    @printf(io," price grid\n")
    for jj in 1:sz.np
        @printf(io,"%16.8f \n",pg[jj])  
    end
    close(io)

    @printf(io," \n")
    @printf(io," transition matrix\n")
    for jj in 1:sz.np
        for ii in 1:sz.np
        @printf(io,"%16.8f",trans[ii,jj]) 
        end
        @printf(io,"\n") 
    end
    close(io)

    y = 1; 
    return y;
end