function interpol(eold::Float64, yold::Float64, aaold::Float64, aold::Float64, dold::Float64,
                  g::NamedTuple, v::Array{Float64,5})
    my_eps = 1.0e-6
    ide = findall(eold .== g.ex)
    idy = findall(yold .== g.y)

    # locate aa (local)
    iaalo = max(findall((aaold - my_eps) .<= g.aa)[1]-1, 1)
    if iaalo == sz.na
        iaahi = sz.na; aafrac = 1.0
    else
        iaahi = iaalo + 1
        aa_hi = g.aa[iaahi]; aa_lo = g.aa[iaalo]
        aafrac = abs(aaold - aa_lo) / abs(aa_hi - aa_lo)
    end

    # locate a (foreign)
    ialo = max(findall((aold - my_eps) .<= g.a)[1]-1, 1)
    if ialo == sz.na
        iahi = sz.na; afrac = 1.0
    else
        iahi = ialo + 1
        a_hi = g.a[iahi]; a_lo = g.a[ialo]
        afrac = abs(aold - a_lo) / abs(a_hi - a_lo)
    end

    # locate d
    iddlo = max(findall((dold - my_eps) .<= g.d)[1]-1, 1)
    if iddlo == sz.nd
        iddhi = sz.nd; dfrac = 1.0
    else
        iddhi = iddlo + 1
        d_hi = g.d[iddhi]; d_lo = g.d[iddlo]
        dfrac = abs(dold - d_lo) / abs(d_hi - d_lo)
    end

    # trilinear over (aa, a, d)
    v000 = v[ide, idy, iaalo, ialo, iddlo]
    v001 = v[ide, idy, iaalo, ialo, iddhi]
    v010 = v[ide, idy, iaalo, iahi, iddlo]
    v011 = v[ide, idy, iaalo, iahi, iddhi]
    v100 = v[ide, idy, iaahi, ialo, iddlo]
    v101 = v[ide, idy, iaahi, ialo, iddhi]
    v110 = v[ide, idy, iaahi, iahi, iddlo]
    v111 = v[ide, idy, iaahi, iahi, iddhi]

    va0 = (1-afrac)*v000 + afrac*v010
    va1 = (1-afrac)*v001 + afrac*v011
    vb0 = (1-afrac)*v100 + afrac*v110
    vb1 = (1-afrac)*v101 + afrac*v111

    v_d0 = (1-aafrac)*va0 + aafrac*vb0
    v_d1 = (1-aafrac)*va1 + aafrac*vb1

    vprime = (1-dfrac)*v_d0 + dfrac*v_d1
    return vprime
end
