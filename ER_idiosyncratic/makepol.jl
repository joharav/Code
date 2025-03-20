function makepol(aint::Array{Int},grid::Vector{Float64});
    enne = size(aint,1);
    enny = size(aint,2);
    enna = size(aint,3);
    ennd = size(aint,4);

    apol = zeros(enne,enny,enna,ennd);
    apg = grid;
    Threads.@threads for id in 1:ennd
        Threads.@threads for ia in 1:enna;
            Threads.@threads for iy in 1:enny;
                Threads.@threads for ie in 1:enne;
                    iia = Int(aint[ie,iy,ia,id]);
                    apol[ie,iy,ia,id] = apg[iia];
                end
            end
        end;
    end;
    return apol::Array{Float64};
end;
