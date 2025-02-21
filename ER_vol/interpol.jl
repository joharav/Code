function interpol(pold::Float64, eold::Float64, aold::Float64, dold::Float64, g::NamedTuple, v::Array{Float64, 4})
    # Extract grids
    pg = g.p
    eg = g.ex
    ag = g.a
    dg = g.d

    # Find the indices for interpolation
    idp = searchsortedfirst(pg, pold)
    ide = searchsortedfirst(eg, eold)
    idalo = searchsortedfirst(ag, aold)
    iddlo = searchsortedfirst(dg, dold)

    # Ensure indices are within bounds
    idp = clamp(idp, 1, length(pg))
    ide = clamp(ide, 1, length(eg))
    idalo = clamp(idalo, 1, length(ag))
    iddlo = clamp(iddlo, 1, length(dg))

    # Calculate the high indices
    idahi = min(idalo + 1, length(ag))
    iddhi = min(iddlo + 1, length(dg))

    # Calculate the fractions for interpolation
    afrac = (aold - ag[idalo]) / (ag[idahi] - ag[idalo])
    dfrac = (dold - dg[iddlo]) / (dg[iddhi] - dg[iddlo])

    # Bilinear interpolation
    v_interp = (1.0 - afrac) * (1.0 - dfrac) * v[idp, ide, idalo, iddlo] +
               afrac * (1.0 - dfrac) * v[idp, ide, idahi, iddlo] +
               (1.0 - afrac) * dfrac * v[idp, ide, idalo, iddhi] +
               afrac * dfrac * v[idp, ide, idahi, iddhi]

    return v_interp
end
