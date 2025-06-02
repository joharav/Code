function makepol(aint::Array{Int},grid::Vector{Float64});
    ennz = size(aint,1);
    enne = size(aint,2);
    enna = size(aint,3);
    ennd = size(aint,4);

    apol = zeros(ennz,enne,enna,ennd);
    apg = grid;
    Threads.@threads for id in 1:ennd
        Threads.@threads for ia in 1:enna;
            Threads.@threads for ie in 1:enne;
                Threads.@threads for iz in 1:ennz;
                    iia = Int(aint[iz,ie,ia,id]);
                    apol[iz,ie,ia,id] = apg[iia];
                end
            end
        end;
    end;
    return apol::Array{Float64};
end;
