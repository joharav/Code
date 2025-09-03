using Plots
function girf_plots(simul_shock::NamedTuple, simul_noshock::NamedTuple)
    # log gap averages across firms
    girf_c  = 100/sz.nFirms * sum(log.(simul_shock.c  ./ simul_noshock.c ), dims=2)
    girf_d  = 100/sz.nFirms * sum(log.(simul_shock.d  ./ simul_noshock.d ), dims=2)
    girf_a  = 100/sz.nFirms * sum(log.(simul_shock.a  ./ simul_noshock.a ), dims=2)
    girf_aa = 100/sz.nFirms * sum(log.(simul_shock.aa ./ simul_noshock.aa), dims=2)

    # effective assets (local units): aa + e*a
    aeff_shock    = simul_shock.aa  .+ simul_shock.ex .* simul_shock.a
    aeff_noshock  = simul_noshock.aa .+ simul_noshock.ex .* simul_noshock.a
    girf_aeff     = 100/sz.nFirms * sum(log.(aeff_shock ./ aeff_noshock), dims=2)

    T_shock = Int(sz.nYears/2)
    rng = T_shock:(T_shock+8)

    # plots
    p1 = plot(girf_c[rng],  title="Consumption response",      xlabel="Quarters", ylabel="% Δ from SS", legend=false)
    savefig(p1, "Output/IRFs/IRF_c.png")

    p2 = plot(girf_d[rng],  title="Durable holdings response", xlabel="Quarters", ylabel="% Δ from SS", legend=false)
    savefig(p2, "Output/IRFs/IRF_d.png")

    p3 = plot(girf_a[rng],  title="Foreign asset response (a)", xlabel="Quarters", ylabel="% Δ from SS", legend=false)
    savefig(p3, "Output/IRFs/IRF_a.png")

    p4 = plot(girf_aa[rng], title="Local asset response (aa)",  xlabel="Quarters", ylabel="% Δ from SS", legend=false)
    savefig(p4, "Output/IRFs/IRF_aa.png")

    p5 = plot(girf_aeff[rng], title="Effective assets response (aa + e·a)", xlabel="Quarters", ylabel="% Δ from SS", legend=false)
    savefig(p5, "Output/IRFs/IRF_a_eff.png")

    return (girf_c, girf_d, girf_a, girf_aa, girf_aeff)
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