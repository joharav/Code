using Printf
using Plots

@inline function safe_log_ratio(x::AbstractArray, y::AbstractArray; eps=1e-10)
    return log.(max.(x, eps)) .- log.(max.(y, eps))
end

function girf_plots(simul_shock::NamedTuple, simul_noshock::NamedTuple; T_shock::Int = Int(sz.nYears ÷ 2))
    outdir = "Output/IRFs"
    isdir(outdir) || mkpath(outdir)

    # log deviations (%), averaged across firms
    girf_c  = 100.0 * vec(mean(safe_log_ratio(simul_shock.c,  simul_noshock.c),  dims=2))
    girf_d  = 100.0 * vec(mean(safe_log_ratio(simul_shock.d,  simul_noshock.d),  dims=2))
    girf_a  = 100.0 * vec(mean(safe_log_ratio(simul_shock.a,  simul_noshock.a),  dims=2))
    girf_aa = 100.0 * vec(mean(safe_log_ratio(simul_shock.aa, simul_noshock.aa), dims=2))

    # effective assets (local units): aa + e*a
    aeff_shock   = simul_shock.aa   .+ simul_shock.ex   .* simul_shock.a
    aeff_noshock = simul_noshock.aa .+ simul_noshock.ex .* simul_noshock.a
    girf_aeff    = 100.0 * vec(mean(safe_log_ratio(aeff_shock, aeff_noshock), dims=2))

    rng = T_shock:min(T_shock + 8, length(girf_c))

    savefig(plot(girf_c[rng],  title="Consumption response", xlabel="Quarters", ylabel="% log dev", legend=false),
            joinpath(outdir, "IRF_c.png"))
    savefig(plot(girf_d[rng],  title="Durables response",    xlabel="Quarters", ylabel="% log dev", legend=false),
            joinpath(outdir, "IRF_d.png"))
    savefig(plot(girf_a[rng],  title="Foreign assets (a)",   xlabel="Quarters", ylabel="% log dev", legend=false),
            joinpath(outdir, "IRF_a.png"))
    savefig(plot(girf_aa[rng], title="Local assets (aa)",    xlabel="Quarters", ylabel="% log dev", legend=false),
            joinpath(outdir, "IRF_aa.png"))
    savefig(plot(girf_aeff[rng], title="Effective assets (aa + e·a)", xlabel="Quarters", ylabel="% log dev", legend=false),
            joinpath(outdir, "IRF_a_eff.png"))

    return (girf_c, girf_d, girf_a, girf_aa, girf_aeff)
end
function compute_cirf(irf::AbstractVector{<:Real}, window_size::Int; name::String="CIRF")
    outdir = "Output/IRFs"
    isdir(outdir) || mkpath(outdir)

    T = length(irf)
    cirf = zeros(Float64, T)
    for t in 1:T
        cirf[t] = sum(@view irf[t:min(T, t + window_size - 1)])
    end

    T_shock = Int(sz.nYears ÷ 2)
    rng = T_shock:min(T_shock + 8, T)
    savefig(plot(cirf[rng], title="Cumulative IRF (window=$window_size)", xlabel="Quarters", ylabel="% log dev", legend=false),
            joinpath(outdir, "$(name)_w$(window_size).png"))

    return cirf
end
