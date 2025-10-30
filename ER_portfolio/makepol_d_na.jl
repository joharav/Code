function makepol_d_na(grid::Vector{Float64}, delta::Float64)
    dpol = zeros(sz.ne, sz.ny, sz.na, sz.na, sz.nd)
    d = grid
    Threads.@threads for id in 1:sz.nd
        Threads.@threads for ia in 1:sz.na
            Threads.@threads for iaa in 1:sz.na
                Threads.@threads for iy in 1:sz.ny
                    Threads.@threads for ie in 1:sz.ne
                        dpol[ie,iy,iaa,ia,id] = (1 - delta) * d[id]
                    end
                end
            end
        end
    end
    return dpol
end
