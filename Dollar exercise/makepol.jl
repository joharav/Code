function makepol(aint::Array{Int},grid::Vector{Float64});
    ennp = size(aint,1);
    enne = size(aint,1);
    ennw = size(aint,1);
    enna = size(aint,2);
    ennd = size(aint,3);

    apol = zeros(ennp,enne,ennw,enna,ennd);
    apg = grid;
    Threads.@threads for id in 1:ennd
        Threads.@threads for ia in 1:enna;
            Threads.@threads for iw in 1:ennw;
                Threads.@threads for ie in 1:enne;
                    Threads.@threads for ip in 1:ennp;
                        iia = Int(aint[ip,ie,iw,ia,id]);
                        apol[ip,ie,iw,ia,id] = apg[iia];
                    end
                end
            end;
        end;
    end;
    return apol::Array{Float64};
end;
