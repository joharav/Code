function makepol_d_na(grid::Vector{Float64}, delta::Float64, chi::Float64)
    dpol = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    d = grid
    Threads.@threads for id in 1:sz.nd
        for ia in 1:sz.na
            for iaa in 1:sz.na
                for iy in 1:sz.ny
                    for ie in 1:sz.ne
                        dpol[ie,iy,iaa,ia,id] = (1 - delta * (1 - chi)) .* d[id];
                    end
                end
            end
        end
    end
    return dpol
end
