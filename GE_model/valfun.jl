function valfun(pea::Vector{Float64})
  # Adjusted value function
  noadjust_result = valfun_noadjust(pea)
  adjust_result = valfun_adjust(pea)

  println("Mean V adjust: ", mean(adjust_result.v))
  println("Mean V noadj: ", mean(noadjust_result.v))

  # Overall value function
  v = max.(adjust_result.v, noadjust_result.v)  # Broadcasting over arrays

  # # Create an indicator matrix (1 if adjusted value is greater, 0 otherwise)
  indicator_matrix = adjust_result.v .> (noadjust_result.v .+ 1e-4)
  adjustment_indicator = vec(indicator_matrix)  # Flattening the matrix to a vector
  num_adjustments = sum(adjustment_indicator.==1)
  total_states = length(adjustment_indicator)
  adjustment_ratio = num_adjustments / total_states
  
  println("Theoretical Adjustment Ratio: $(round(adjustment_ratio * 100, digits=2))%")

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
  outtuple = (v = v, gidx = gidx, pol = pol, g = g, e = e, adjust_result=adjust_result, noadjust_result=noadjust_result,adjustment_indicator=indicator_matrix) 

  # Optional debug or visualization
  if settings.verbose && adjust_result.e == 0 && noadjust_result.e == 0
      printstuff(outtuple)
      plotstuff(outtuple.v,gidx_a,gidx_d,outtuple.pol.c,g) 
  end

  return outtuple
end
