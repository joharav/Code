function valfun(pea::Vector{Float64})
  # Adjusted value function
  adjust_result = valfun_adjust(pea)
  noadjust_result = valfun_noadjust(pea)

  # Overall value function
  v = max.(adjust_result.v, noadjust_result.v)  # Broadcasting over arrays

  # # Create an indicator matrix (1 if adjusted value is greater, 0 otherwise)
  indicator_matrix = adjust_result.v .> noadjust_result.v
  num_adjustments = sum(adjust_result.v .> noadjust_result.v)
  total_states = length(adjust_result.v)
  adjustment_ratio = num_adjustments / total_states
  
  println("Total Adjustments: ", num_adjustments)
  println("Total States: ", total_states)
  println("Adjustment Ratio: ", adjustment_ratio)

  
  # # Conditional assignments based on the indicator_matrix
  gidx_a = ifelse.(indicator_matrix, adjust_result.gidx.a, noadjust_result.gidx.a)
  gidx_d = ifelse.(indicator_matrix, adjust_result.gidx.d, noadjust_result.gidx.d)
  gidx = dtp.Ipol(gidx_a, gidx_d)  # Reconstructing the Ipol object

  pol_a = ifelse.(indicator_matrix, adjust_result.pol.a, noadjust_result.pol.a)
  pol_d = ifelse.(indicator_matrix, adjust_result.pol.d, noadjust_result.pol.d)
  pol = dtp.Pol(pol_a, pol_d)  # Reconstructing the Pol object

   #mew=ifelse.(indicator_matrix, adjust_result.mew, noadjust_result.mew)

  g = adjust_result.g
  e = adjust_result.e

  # Constructing the output tuple
  outtuple = (v = v, gidx = gidx, pol = pol, g = g, e = e, adjust_result=adjust_result, noadjust_result=noadjust_result) #, mew = mew

  # Optional debug or visualization
  if settings.verbose && adjust_result.e == 0 && noadjust_result.e == 0
      printstuff(outtuple)
      plotstuff(outtuple.v,gidx_a,gidx_d,g) #, mew
  end

  return outtuple
end
