function bootstrap_girf(simul_shock, simul_noshock, B=1000, alpha=0.05)
    girf_bootstrap = zeros(sz.nYears, B) 

    for b in 1:B
        resample_idx = rand(1:sz.nFirms, sz.nFirms)  # Sample firms with replacement
        girf_resample = 100/sz.nFirms * sum(simul_shock[:, resample_idx] .- simul_noshock[:, resample_idx], dims=2)
        girf_bootstrap[:, b] = girf_resample[:]
    end

    lower_bound = mapslices(x -> quantile(x, alpha/2), girf_bootstrap, dims=2)
    upper_bound = mapslices(x -> quantile(x, 1 - alpha/2), girf_bootstrap, dims=2)
    mean_girf   = mean(girf_bootstrap, dims=2)

    outuple=(mean=mean_girf, lower=lower_bound, upper=upper_bound)

    return outuple
end

function plot_irf_with_ci(outuple::NamedTuple, title_str::String, filename::String)
    irf_mean=outuple.mean
    irf_lower=outuple.lower
    irf_upper=outuple.upper

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
    girf_c = bootstrap_girf(simul_shock.c, simul_noshock.c)
    girf_d = bootstrap_girf(simul_shock.d, simul_noshock.d)
    girf_a = bootstrap_girf(simul_shock.a, simul_noshock.a)

    # Plot each IRF
    plot_irf_with_ci(girf_c, "IRF - c", "Output/IRFs/IRF_c.png")
    plot_irf_with_ci(girf_d, "IRF - d", "Output/IRFs/IRF_d.png")
    plot_irf_with_ci(girf_a, "IRF - a", "Output/IRFs/IRF_a.png")


    cirf_c          =compute_cirf(girf_c.mean, 8,"c");
    cirf_d          =compute_cirf(girf_d.mean, 8,"d");
    cirf_a          =compute_cirf(girf_a.mean, 8,"a");


    outuple = (girf_c, girf_d, girf_a, cirf_c, cirf_d, cirf_a)

    return outuple
end

# Function to compute cumulative IRF over rolling windows
function compute_cirf(irf_results::Vector{Float64}, window_size::Int64,x::String)
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
