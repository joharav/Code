function printstuff(answ::NamedTuple)

    vnew = answ.v

    pol = answ.pol
    apol = pol.a
    dpol = pol.d
    grids = answ.g

    filename = "Exchange_Rate/Output/v.txt"
    io = open(filename, "w")
    
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    @printf(io, "%16.8f", vnew[ip, ie, ia, id])
                    if ip < sz.np
                        @printf(io, " ")  # Space between p values
                    end
                end
                @printf(io, "\n")  # New line for each e within a w block
            end
            @printf(io, "\n")  # New line for each a within a d block
        end
        @printf(io, "\n")  # Extra newline to separate blocks of d for readability
    end
    
    close(io)


    filename = "Exchange_Rate/Output/a.txt"
    io = open(filename, "w")
    
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    @printf(io, "%16.8f", apol[ip, ie, ia, id])
                    if ip < sz.np
                        @printf(io, " ")  # Space between p values
                    end
                end
                @printf(io, "\n")  # New line for each e within a w block
            end
            @printf(io, "\n")  # New line for each w within an a block
        end
        @printf(io, "\n")  # Extra newline to separate blocks of d for readability
    end
    
    close(io)


    filename = "Exchange_Rate/Output/d.txt"
    io = open(filename, "w")
    
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for ie in 1:sz.ne
                Threads.@threads for ip in 1:sz.np
                    @printf(io, "%16.8f", dpol[ip, ie, ia, id])
                    if ip < sz.np
                        @printf(io, " ")  # Space between p values
                    end
                end
                @printf(io, "\n")  # New line for each e within a w block
            end
            @printf(io, "\n")  # New line for each w within an a block
        end
        @printf(io, "\n")  # Extra newline to separate blocks of d for readability
    end
    
    close(io)

    pg = grids.p;
    eg = grids.ex;
    ag = grids.a;
    apg = grids.ap;
    dg = grids.d;
    dpg = grids.dp;
    trans = grids.t; 
    filname = "Exchange_Rate/Output/statespace.txt"
    io = open(filname,"w")
    @printf(io," asset state grid\n")
    Threads.@threads for jj in 1:sz.na
        @printf(io,"%16.8f \n",ag[jj])  
    end
    @printf(io," \n")
    @printf(io," asset policy grid\n")
    Threads.@threads for jj in 1:sz.npa
        @printf(io,"%16.8f \n",apg[jj])  
    end
    @printf(io," \n")
    @printf(io," durable state grid\n")
    Threads.@threads for jj in 1:sz.nd
        @printf(io,"%16.8f \n",dg[jj])  
    end
    @printf(io," \n")
    @printf(io," durable policy grid\n")
    Threads.@threads for jj in 1:sz.npd
        @printf(io,"%16.8f \n",dpg[jj])  
    end
    @printf(io," price grid\n")
    Threads.@threads for jj in 1:sz.np
        @printf(io,"%16.8f \n",pg[jj])  
    end
    @printf(io," \n")
    @printf(io," exchange rate grid\n")
    Threads.@threads for jj in 1:sz.ne
        @printf(io,"%16.8f \n",eg[jj])  
    end
    @printf(io," \n")
    @printf(io," transition matrix\n")
    Threads.@threads for jj in 1:(sz.np * sz.ne)
        Threads.@threads for ii in 1:(sz.np * sz.ne)
            @printf(io,"%16.8f",trans[ii,jj]) 
        end
        @printf(io,"\n") 
    end
    close(io)
    
    y = 1; 
    return y;
end