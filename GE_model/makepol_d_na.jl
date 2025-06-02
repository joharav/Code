function makepol_d_na(grid::Vector{Float64});
    delta       = pea[2]        # Depreciation rate for durables
    chi         = pea[9]        # Required maintenance

    dpol = zeros(sz.nz,sz.ne,sz.na,sz.nd);
    d = grid;

    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na;
            Threads.@threads for ie in 1:sz.ne;
                Threads.@threads for iz in 1:sz.nz;
                   dpol[iz,ie,ia,id] = (1 - delta * (1 - chi)) .* d[id];
                end
            end
        end;
    end;
    return dpol::Array{Float64};
end;
