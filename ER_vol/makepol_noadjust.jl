function makepol_noadjust(gidx::dtp.Ipol, grids::NamedTuple)
    delta = pea[2]    # Depreciation rate for durables
    chi = pea[11]     # Required maintenance 

    pol_a = grids.ap[gidx.a]
    pol_d = (1 - delta * (1 - chi)) * grids.d[gidx.d]
    return dtp.Pol(pol_a, pol_d)
end