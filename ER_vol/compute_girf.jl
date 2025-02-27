function(simul_shock::NamedTuple,simul_noshock::NamedTuple)
    girf_v          = 100/sz.nFirms * sum(log.(simul_shock.allv_shock ./ simul_noshock.allv_noshock), dims=2)
    girf_c          = 100/sz.nFirms * sum(log.(simul_shock.allc_shock ./ simul_noshock.allc_noshock), dims=2)
    girf_d          = 100/sz.nFirms * sum(log.(simul_shock.alld_shock ./ simul_noshock.alld_noshock), dims=2)
    girf_a          = 100/sz.nFirms * sum(log.(simul_shock.alla_shock ./ simul_noshock.alla_noshock), dims=2)
    girf_dadjust    = 100/sz.nFirms * sum(log.(simul_shock.alld_adjust_shock ./ simul_noshock.alld_adjust_noshock), dims=2)
    girf_adjust=(100/sz.nFirms)*sum(simul_shock.alld_adjust_shock.-simul_noshock.alld_adjust_noshock, dims=2)

    #Plot IRFs
    plot1=plot(girf_v, title="IRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot1, "Output/IRFs/IRF_v.png")

    plot2=plot(girf_c, title="IRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot2, "Output/IRFs/IRF_c.png")

    plot3=plot(girf_d, title="IRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot3, "Output/IRFs/IRF_d.png")

    plot4=plot(girf_a, title="IRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot4, "Output/IRFs/IRF_a.png")

    plot5=plot(girf_dadjust, title="IRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot5, "Output/IRFs/IRF_dadjust.png")

    plot6=plot(girf_adjust, title="IRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot6, "Output/IRFs/IRF_adjust.png")

    outuple=(girf_v = girf_v, girf_c = girf_c, girf_d = girf_d, girf_a = girf_a, girf_dadjust = girf_dadjust, girf_adjust = girf_adjust)

    return outuple
end


# Function to compute cumulative IRF over rolling windows
function compute_cirf(irf_results::Vector{Float64}, window_size::Int,x::String)
    param_str = string(x, "_", @sprintf("%.4f", window_size))

    irf_periods = length(irf_results)
    cirf_vector = zeros(irf_periods)

    for t in 1:irf_periods
        start_period = t
        end_period = min(t + window_size - 1, irf_periods)  # Ensure we don't go out of bounds
        cirf_vector[t] = sum(irf_results[start_period:end_period])
    end

    plot=plot(cirf_vector, title="Cumulative CIRF", xlabel="Quarters", ylabel="% Change from SS")
    savefig(plot, "Output/IRFs/CIRF_($param_str).png")


    return cirf_vector
end