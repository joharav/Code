function makepol(aint::Array{Int},grid::Vector{Float64});
    enne = size(aint,1);
    enna = size(aint,2);
    ennd = size(aint,3);

    apol = zeros(enne,enna,ennd);
    apg = grid;
    Threads.@threads for id in 1:ennd
        Threads.@threads for ia in 1:enna;
            Threads.@threads for ie in 1:enne;
                iia = Int(aint[ie,ia,id]);
                apol[ie,ia,id] = apg[iia];
            end
        end;
    end;
    return apol::Array{Float64};
end;
