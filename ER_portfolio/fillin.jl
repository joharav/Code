# obj: (ne, ny, naa, na, nd)  ->  objlong: (ne, ny, npa, npa, npd)
function fillin(obj::Array{Float64,5}, g::NamedTuple)
  objlong = zeros(sz.ne, sz.ny, sz.npa, sz.npa, sz.npd)

  # Thread over all (iiaa, iia, iid) combos via a single linear index
  Threads.@threads for J in 1:(sz.npa * sz.npa * sz.npd)
      iiaa = 1 + (J - 1) % sz.npa
      tmp  = (J - 1) ÷ sz.npa
      iia  = 1 + tmp % sz.npa
      iid  = 1 + tmp ÷ sz.npa

      # ----- aa-weights (state aa -> policy aap) -----
      aa_down = Int(floor((sz.na - 1.0) * (iiaa - 1.0) / (sz.npa - 1)) + 1)
      aa_up   = (iiaa == sz.npa) ? aa_down : aa_down + 1
      den_aa  = g.aa[aa_up] - g.aa[aa_down]
      waa = (iiaa == sz.npa || den_aa == 0.0) ? 1.0 : (g.aap[iiaa] - g.aa[aa_down]) / den_aa

      # ----- a-weights (state a -> policy ap) -----
      a_down = Int(floor((sz.na - 1.0) * (iia - 1.0) / (sz.npa - 1)) + 1)
      a_up   = (iia == sz.npa) ? a_down : a_down + 1
      den_a  = g.a[a_up] - g.a[a_down]
      wa = (iia == sz.npa || den_a == 0.0) ? 1.0 : (g.ap[iia] - g.a[a_down]) / den_a

      # ----- d-weights (state d -> policy dp) -----
      d_down = Int(floor((sz.nd - 1.0) * (iid - 1.0) / (sz.npd - 1)) + 1)
      d_up   = (iid == sz.npd) ? d_down : d_down + 1
      den_d  = g.d[d_up] - g.d[d_down]
      wd = (iid == sz.npd || den_d == 0.0) ? 1.0 : (g.dp[iid] - g.d[d_down]) / den_d

      # ----- 3D trilinear over (aa, a, d) for each (ie,iy) -----
      for iy in 1:sz.ny
          for ie in 1:sz.ne
              v000 = obj[ie, iy, aa_down, a_down, d_down]
              v001 = obj[ie, iy, aa_down, a_down, d_up  ]
              v010 = obj[ie, iy, aa_down, a_up,   d_down]
              v011 = obj[ie, iy, aa_down, a_up,   d_up  ]
              v100 = obj[ie, iy, aa_up,   a_down, d_down]
              v101 = obj[ie, iy, aa_up,   a_down, d_up  ]
              v110 = obj[ie, iy, aa_up,   a_up,   d_down]
              v111 = obj[ie, iy, aa_up,   a_up,   d_up  ]

              va0 = (1 - wa) * v000 + wa * v010
              va1 = (1 - wa) * v001 + wa * v011
              vb0 = (1 - wa) * v100 + wa * v110
              vb1 = (1 - wa) * v101 + wa * v111

              v_d0 = (1 - waa) * va0 + waa * vb0
              v_d1 = (1 - waa) * va1 + waa * vb1

              objlong[ie, iy, iiaa, iia, iid] = (1 - wd) * v_d0 + wd * v_d1
          end
      end
  end
  return objlong
end
