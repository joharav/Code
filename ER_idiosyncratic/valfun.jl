function valfun(pea::Vector{Float64})
  # Adjusted value function
  noadjust_result = valfun_noadjust(pea)
  adjust_result = valfun_adjust(pea)

  # Overall value function
  v = max.(adjust_result.v, noadjust_result.v)  # Broadcasting over arrays

  # # Create an indicator matrix (1 if adjusted value is greater, 0 otherwise)
  indicator_matrix = adjust_result.v .> noadjust_result.v
  adjust_flag = Float64.(adjust_result.v .> noadjust_result.v)

  #print how many firms are adjusted
  if settings.verbose
      println("Number of firms with adjusted value function: ", 100*sum(indicator_matrix)/length(indicator_matrix))
  end

  # # Conditional assignments based on the indicator_matrix
  gidx_a = ifelse.(indicator_matrix, adjust_result.gidx.a, noadjust_result.gidx.a)
  gidx_d = ifelse.(indicator_matrix, adjust_result.gidx.d, noadjust_result.gidx.d)
  gidx = dtp.Ipol(gidx_a, gidx_d)  # Reconstructing the Ipol object

  pol_a = ifelse.(indicator_matrix, adjust_result.pol.a, noadjust_result.pol.a)
  pol_d = ifelse.(indicator_matrix, adjust_result.pol.d, noadjust_result.pol.d)
  pol_c = ifelse.(indicator_matrix, adjust_result.pol.c, noadjust_result.pol.c)

  pol = dtp.Pol(pol_a, pol_d, pol_c)  # Reconstructing the Pol object

  g = adjust_result.g
  e = adjust_result.e

  # Constructing the output tuple
  outtuple = (v = v, gidx = gidx, pol = pol, g = g, e = e, adjust_result=adjust_result, noadjust_result=noadjust_result, adjust_flag = adjust_flag) 

  # Optional debug or visualization
  if settings.verbose && adjust_result.e == 0 && noadjust_result.e == 0
      printstuff(outtuple)
      plotstuff(outtuple.v,gidx_a,gidx_d,outtuple.pol.c,g) 
  end

  return outtuple
end
