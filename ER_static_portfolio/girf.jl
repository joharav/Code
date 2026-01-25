@inline function _hasfield(nt::NamedTuple, s::Symbol)
        return s in propertynames(nt)
    end
    
    @inline function _derive_a_aa_from_wse(w, s, e; eps=1e-12)
        aa = (1 .- s) .* w
        a  = (s .* w) ./ max.(e, eps)
        return a, aa
    end
    
    @inline function safe_log_ratio(x::AbstractArray, y::AbstractArray; eps=1e-10)
        return log.(max.(x, eps)) .- log.(max.(y, eps))
    end
    
    function girf_plots(sim_shock::NamedTuple, sim_base::NamedTuple; T_shock::Int = Int(sz.nYears ÷ 2))
        outdir = "Output/IRFs"
        isdir(outdir) || mkpath(outdir)
    
        # consumption & durables always exist
        girf_c = 100.0 * vec(mean(safe_log_ratio(sim_shock.c, sim_base.c), dims=2))
        girf_d = 100.0 * vec(mean(safe_log_ratio(sim_shock.d, sim_base.d), dims=2))
    
        # foreign/local assets: either use stored a/aa or derive from (w,s,e)
        a_sh, aa_sh = _hasfield(sim_shock, :a) && _hasfield(sim_shock, :aa) ?
            (sim_shock.a, sim_shock.aa) :
            _derive_a_aa_from_wse(sim_shock.w, sim_shock.s, sim_shock.ex)
    
        a_b, aa_b = _hasfield(sim_base, :a) && _hasfield(sim_base, :aa) ?
            (sim_base.a, sim_base.aa) :
            _derive_a_aa_from_wse(sim_base.w, sim_base.s, sim_base.ex)
    
        girf_a  = 100.0 * vec(mean(safe_log_ratio(a_sh,  a_b),  dims=2))
        girf_aa = 100.0 * vec(mean(safe_log_ratio(aa_sh, aa_b), dims=2))
    
        # effective assets (local): aa + e*a
        aeff_sh = aa_sh .+ sim_shock.ex .* a_sh
        aeff_b  = aa_b  .+ sim_base.ex  .* a_b
        girf_aeff = 100.0 * vec(mean(safe_log_ratio(aeff_sh, aeff_b), dims=2))
    
        rng = T_shock:min(T_shock + 8, length(girf_c))
    
        savefig(plot(girf_c[rng], xlabel="Quarters", ylabel="% log dev", title="Consumption response", legend=false),
                joinpath(outdir, "IRF_c.pdf"))
        savefig(plot(girf_d[rng], xlabel="Quarters", ylabel="% log dev", title="Durables response", legend=false),
                joinpath(outdir, "IRF_d.pdf"))
        savefig(plot(girf_a[rng], xlabel="Quarters", ylabel="% log dev", title="Foreign assets (a)", legend=false),
                joinpath(outdir, "IRF_a.pdf"))
        savefig(plot(girf_aa[rng], xlabel="Quarters", ylabel="% log dev", title="Local assets (aa)", legend=false),
                joinpath(outdir, "IRF_aa.pdf"))
        savefig(plot(girf_aeff[rng], xlabel="Quarters", ylabel="% log dev", title="Effective assets (aa + e·a)", legend=false),
                joinpath(outdir, "IRF_a_eff.pdf"))
    
        return (girf_c, girf_d, girf_a, girf_aa, girf_aeff)
    end
    
    function compute_cirf(irf::AbstractVector{<:Real}, window_size::Int; name::String="CIRF", T_shock::Int=Int(sz.nYears ÷ 2))
        outdir = "Output/IRFs"
        isdir(outdir) || mkpath(outdir)
    
        T = length(irf)
        cirf = zeros(Float64, T)
        for t in 1:T
            cirf[t] = sum(@view irf[t:min(T, t + window_size - 1)])
        end
    
        rng = T_shock:min(T_shock + 8, T)
        savefig(plot(cirf[rng], title="Cumulative IRF (window=$window_size)", xlabel="Quarters", ylabel="% log dev", legend=false),
                joinpath(outdir, "$(name)_w$(window_size).pdf"))
    
        return cirf
    end
    
