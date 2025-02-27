function bootstrap_girf(simul_shock, simul_noshock, B=1000, alpha=0.05)
    nFirms = size(simul_shock.d, 2)
    nYears = size(simul_shock.d, 1)
    girf_bootstrap = zeros(nYears, B) 

    for b in 1:B
        resample_idx = rand(1:nFirms, nFirms)  # Sample firms with replacement
        girf_resample = 100/nFirms * sum(log.(simul_shock.d[:, resample_idx] ./ simul_noshock.d[:, resample_idx]), dims=2)
        girf_bootstrap[:, b] = girf_resample[:]
    end

    lower_bound = quantile(girf_bootstrap, alpha/2, dims=2)
    upper_bound = quantile(girf_bootstrap, 1 - alpha/2, dims=2)
    mean_girf   = mean(girf_bootstrap, dims=2)

    return (mean=mean_girf, lower=lower_bound, upper=upper_bound)
end

function plot_irf_with_ci(irf_mean::Vector{Float64}, irf_lower::Vector{Float64}, 
    irf_upper::Vector{Float64}, title_str::String, filename::String)
    x_axis = 1:length(irf_mean)  # X-axis representing time periods (quarters)

    plot1 = plot(x_axis, irf_mean, lw=2, label="Mean IRF", title=title_str, xlabel="Quarters", ylabel="% Change from SS")
    plot!(x_axis, irf_lower, lw=1, linestyle=:dash, label="Lower Bound", color=:gray)
    plot!(x_axis, irf_upper, lw=1, linestyle=:dash, label="Upper Bound", color=:gray)
    
    # Shaded confidence interval
    plot!(x_axis, irf_upper, fill_between=irf_lower, color=:gray, fillalpha=0.3, label="Confidence Interval")

    savefig(plot1, filename)
end

function compute_and_plot_irfs(simul_shock, simul_noshock)
    # Compute GIRFs with bootstrap CIs
    girf_v = bootstrap_girf(simul_shock.allv_shock, simul_noshock.allv_noshock)
    girf_c = bootstrap_girf(simul_shock.allc_shock, simul_noshock.allc_noshock)
    girf_d = bootstrap_girf(simul_shock.alld_shock, simul_noshock.alld_noshock)
    girf_a = bootstrap_girf(simul_shock.alla_shock, simul_noshock.alla_noshock)
    girf_dadjust = bootstrap_girf(simul_shock.alld_adjust_shock, simul_noshock.alld_adjust_noshock)

    # Plot each IRF
    plot_irf_with_ci(girf_v.mean, girf_v.lower, girf_v.upper, "IRF - v", "Output/IRFs/IRF_v.png")
    plot_irf_with_ci(girf_c.mean, girf_c.lower, girf_c.upper, "IRF - c", "Output/IRFs/IRF_c.png")
    plot_irf_with_ci(girf_d.mean, girf_d.lower, girf_d.upper, "IRF - d", "Output/IRFs/IRF_d.png")
    plot_irf_with_ci(girf_a.mean, girf_a.lower, girf_a.upper, "IRF - a", "Output/IRFs/IRF_a.png")
    plot_irf_with_ci(girf_dadjust.mean, girf_dadjust.lower, girf_dadjust.upper, "IRF - dadjust", "Output/IRFs/IRF_dadjust.png")

    return (girf_v, girf_c, girf_d, girf_a, girf_dadjust)
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
