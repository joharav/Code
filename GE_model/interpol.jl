function interpol(zold::Float64, eold::Float64, aold::Float64, dold::Float64, g::NamedTuple, v::Array{Float64})
    my_eps = 1.0e-6
    # Exact match for e, as before
    ide = findall(eold .== g.ex);  # easy because it is an exact match
    idy = findall(zold .== g.zz);  # easy because it is an exact match

    # Interpolation in a dimension
    idalo = findall((aold - my_eps) .<= g.a);
    idalo = max(idalo[1]-1, 1); # Ensure g.k is sorted in ascending order
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
    iddlo = findall((dold - my_eps) .<= g.d);#println(idplo)
    iddlo = max(iddlo[1]-1, 1); # Similar logic for p as for k
    if iddlo == sz.nd
        iddhi = sz.nd
        dfrac = 1.0
    else
        iddhi = iddlo + 1
        dhi = g.d[Int(iddhi)]
        dlo = g.d[Int(iddlo)]
        dfrac = abs(dold - dlo) / abs(dhi - dlo)
    end
    
    # Interpolating in both k and p dimensions
    v11 = v[ide, idy, idalo, iddlo]
    v12 = v[ide, idy, idalo, iddhi]
    v21 = v[ide, idy, idahi, iddlo]
    v22 = v[ide, idy, idahi, iddhi]
   
    # Bilinear interpolation formula
    vprime = afrac * (dfrac * v22 + (1 - dfrac) * v12) +
             (1 - afrac) * (dfrac * v21 + (1 - dfrac) * v11)
            
            
   
    return vprime
end
