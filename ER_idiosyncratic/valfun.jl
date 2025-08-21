using StatsBase
function valfun(pea::Vector{Float64}; λ::Float64 = 0.0)
  # Adjusted value function
  noadjust_result = valfun_noadjust(pea; λ=λ)
  if settings.verbose
    println("Non-adjust: mean(v) = ", mean(noadjust_result.v))
    println("Non-adjust: min(v) = ", minimum(noadjust_result.v))
    println("Non-adjust: max(v) = ", maximum(noadjust_result.v))
  end
 
  

  adjust_result = valfun_adjust(pea; λ=λ)
  if settings.verbose
    println("Adjust: mean(v) = ", mean(adjust_result.v))
    println("Adjust: min(v) = ", minimum(adjust_result.v))
    println("Adjust: max(v) = ", maximum(adjust_result.v))
  end
  # Overall value function
  v = max.(adjust_result.v, noadjust_result.v)  # Broadcasting over arrays

  val_diff = adjust_result.v .- noadjust_result.v

  if settings.verbose
    println("Mean diff: ", mean(val_diff), " | Min: ", minimum(val_diff), " | Max: ", maximum(val_diff))
    println("Share where adjust gives lower value: ", sum(val_diff .< 0) / length(val_diff))
  end

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

  if settings.verbose
    println("Mean v adjust: ", mean(adjust_result.v))
    println("Mean v noadjust: ", mean(noadjust_result.v))
    println("Mean a adjust: ", mean(adjust_result.pol.a))
    println("Mean a noadjust: ", mean(noadjust_result.pol.a))
    println("Mean d adjust: ", mean(adjust_result.pol.d))
    println("Mean d noadjust: ", mean(noadjust_result.pol.d))
  end


  # Constructing the output tuple
  outtuple = (v = v, gidx = gidx, pol = pol, g = g, e = e, adjust_result=adjust_result, noadjust_result=noadjust_result, adjust_flag = adjust_flag,adjustment_indicator=indicator_matrix) 

  # Optional debug or visualization
  if settings.verbose && adjust_result.e == 0 && noadjust_result.e == 0
      printstuff(outtuple)
      plotstuff(outtuple.v,gidx_a,gidx_d,outtuple.pol.c,g) 
  end

  return outtuple
end
