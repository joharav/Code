function makepol_noadjust(gidx::dtp.Ipol, grids::NamedTuple)
    pol_a = grids.ap[gidx.a]
    pol_d = (1 - grids.delta * (1 - grids.chi)) * grids.d[gidx.d]
    return dtp.Pol(pol_a, pol_d)
end