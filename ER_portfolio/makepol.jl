# === new 5D version for states (ie, iy, iaa, ia, id) ===
function makepol(aint::Array{Int,5}, grid::Vector{Float64})
    enne, enny, ennaa, enna, ennd = size(aint)
    apol = zeros(enne, enny, ennaa, enna, ennd)
    apg  = grid
    Threads.@threads for id in 1:ennd
        for ia in 1:enna
            for iaa in 1:ennaa
                for iy in 1:enny
                    for ie in 1:enne
                        iix = Int(aint[ie,iy,iaa,ia,id])
                        apol[ie,iy,iaa,ia,id] = apg[iix]
                    end
                end
            end
        end
    end
    return apol::Array{Float64,5}
end
