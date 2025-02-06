function fillin(obj::Array{Float64}, g::NamedTuple)

  # Assuming obj is now a three-dimensional array,
  # and we're interpolating along the last two dimensions (k and p).
  objlong = zeros(sz.np, sz.npa, sz.npd)

  Threads.@threads for ia in 1:sz.npa
      for id in 1:sz.npd
          # Interpolation in the k dimension
          adown = Int(floor((sz.na-1.0) * (ia-1.0) / (sz.npa-1)) + 1)
          aup = ia == sz.npa ? adown : adown + 1 #slick way of doing if else endif
          afrac = ia == sz.npa ? 1.0 : (g.ap[ia] - g.a[adown]) / (g.a[aup] - g.a[adown])

          # Interpolation in the p dimension
          ddown = Int(floor((sz.nd-1.0) * (id-1.0) / (sz.npd-1)) + 1)
          dup = id == sz.npd ? ddown : ddown + 1
          dfrac = id == sz.npd ? 1.0 : (g.dp[id] - g.d[ddown]) / (g.d[dup] - g.d[ddown])

          # Two-dimensional interpolation
          for ip in 1:sz.np
              objlong[ip, ia, id] = afrac * (dfrac * obj[ip, aup, dup] + (1 - dfrac) * obj[ip, aup, ddown]) +
                                    (1 - afrac) * (dfrac * obj[ip, adown, dup] + (1 - dfrac) * obj[ip, adown, ddown])
          end
      end
  end

  return objlong::Array{Float64}
end
