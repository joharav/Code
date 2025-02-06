function makepol(aint::Array{Int},grid::Vector{Float64});
    ennp = size(aint,1);
    enna = size(aint,2);
    ennd = size(aint,3);

    apol = zeros(ennp,enna,ennd);
    apg = grid;
    Threads.@threads for id in 1:ennd
        for ia in 1:enna;
            for ip in 1:ennp;
               iia = Int(aint[ip,ia,id]);
               apol[ip,ia,id] = apg[iia];
            end;
        end;
    end;
    return apol::Array{Float64};
end;
