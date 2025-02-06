function interpol(pold::Float64, eold::Float64, wold::Float64, aold::Float64, dold::Float64, g::NamedTuple, v::Array{Float64})
    my_eps = 1.0e-6

    # Exact match for p, e, and w
    idp = findfirst(pold .== g.p)  # easy because it is an exact match
    ide = findfirst(eold .== g.e)  # easy because it is an exact match
    idw = findfirst(wold .== g.w)  # easy because it is an exact match

    # Interpolation in a dimension
    idalo = findfirst((aold - my_eps) .<= g.a)
    idalo = max(idalo-1, 1)  # Ensure g.a is sorted in ascending order
    if idalo == sz.na
        idahi = sz.na
        afrac = 1.0
    else
        idahi = idalo + 1
        ahi = g.a[Int(idahi)]
        alo = g.a[Int(idalo)]
        afrac = abs(aold - alo) / abs(ahi - alo)
    end

    # Interpolation in d dimension
    iddlo = findfirst((dold - my_eps) .<= g.d)
    iddlo = max(iddlo-1, 1)  # Similar logic for d as for a
    if iddlo == sz.nd
        iddhi = sz.nd
        dfrac = 1.0
    else
        iddhi = iddlo + 1
        dhi = g.d[Int(iddhi)]
        dlo = g.d[Int(iddlo)]
        dfrac = abs(dold - dlo) / abs(dhi - dlo)
    end

    # Bilinear interpolation
    v_interp = (1.0 - afrac) * (1.0 - dfrac) * v[idp, ide, idw, idalo, iddlo] +
               afrac * (1.0 - dfrac) * v[idp, ide, idw, idahi, iddlo] +
               (1.0 - afrac) * dfrac * v[idp, ide, idw, idalo, iddhi] +
               afrac * dfrac * v[idp, ide, idw, idahi, iddhi]

    return v_interp
end
