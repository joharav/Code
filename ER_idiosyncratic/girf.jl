using Plots
function girf_plots(simul_shock::NamedTuple,simul_noshock::NamedTuple)
     girf_c          = 100/sz.nFirms * sum(log.(simul_shock.c ./ simul_noshock.c), dims=2)
     girf_d          = 100/sz.nFirms * sum(log.(simul_shock.d ./ simul_noshock.d), dims=2)
     girf_a          = 100/sz.nFirms * sum(log.(simul_shock.a ./ simul_noshock.a), dims=2)

     T_shock=Int(sz.nYears/2)
     #Plot IRFs
     plot1=plot(girf_c[T_shock:T_shock+8], title="Consumption response", xlabel="Quarters", ylabel="% Change from SS", legend=false)
     savefig(plot1, "Output/IRFs/IRF_c.png")
 
     plot2=plot(girf_d[T_shock:T_shock+8], title="Durable holdings response", xlabel="Quarters", ylabel="% Change from SS", legend=false)
     savefig(plot2, "Output/IRFs/IRF_d.png")
     plot3=plot(girf_a[T_shock:T_shock+8], title="Asset holdings response", xlabel="Quarters", ylabel="% Change from SS", legend=false)
     savefig(plot3, "Output/IRFs/IRF_a.png")
 
     outuple=( girf_c, girf_d, girf_a)
 
     return outuple
 end
 
 
# Function to compute cumulative IRF over rolling windows
function compute_cirf(irf_results::Vector{Float64}, window_size::Int,x::String)
    param_str = string(x, "_", @sprintf("%.4f", window_size))

    irf_periods = length(irf_results)
    cirf_vector = zeros(irf_periods)
    T_shock=Int(sz.nYears/2)

    for t in 1:irf_periods
        start_period = t
        end_period = min(t + window_size - 1, irf_periods)  # Ensure we don't go out of bounds
        cirf_vector[t] = sum(irf_results[start_period:end_period])
    end

    pp1=plot(cirf_vector[T_shock:T_shock+8], title="Cumulative CIRF", xlabel="Quarters", ylabel="% Change from SS", legend=false)
    savefig(pp1, "Output/IRFs/CIRF_($param_str).png")


    return cirf_vector
end