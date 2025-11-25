using Statistics

# ---------- helpers ----------
_fin(v) = v[isfinite.(v)]
@inline _safe_den(x) = max(x, 1e-12)


# Running time-since-last-adjustment (in periods) for each (t,i) in r0
# adj: T×N Bool/Int in the selected window r0
function running_spells(adj::AbstractMatrix)
    T, N = size(adj)
    out = zeros(Float64, T, N)
    @inbounds for j in 1:N
        s = 0.0
        for t in 1:T
            if adj[t,j] > 0
                s = 0.0            # reset on adjustment
            else
                s += 1.0
            end
            out[t,j] = s
        end
    end
    return out  # same size as adj
end





function makemoments(simdata::NamedTuple, pea::Vector{Float64};
    per_year::Int = 4, shock::Bool=false)

# --- parameters ---
w  = pea[8]    # earnings scale
pd = pea[10]   # durable price (USD)

# --- post-burn-in window ---
r0 = (sz.burnin-2):sz.nYears
r1 = (sz.burnin-3):(sz.nYears-1)

# slice panel
a      = simdata.a[r0, :]
aa     = simdata.aa[r0, :]
d      = simdata.d[r0, :]
ex     = simdata.ex[r0, :]
y      = simdata.y[r0, :]
adj    = simdata.adjust_indicator[r0, :]

# --- value/denominators in local units, panel-wide ---
Vd    = pd .* ex .* d                         # durable value
a_loc = aa
a_fx  = ex .* a                               # USD assets valued in local
a_eff = a_loc .+ a_fx                         # total liquid in local

# --- build panel vectors (flatten time×households) ---
d_value_vec     = vec(Vd)
a_eff_vec       = vec(a_eff)
a_fx_vec        = vec(a_fx)

# shares/ratios at (t,i)
d_wealth_ratio  = d_value_vec ./ max.(d_value_vec .+ a_eff_vec, 1e-12)
usd_sh          = a_fx_vec    ./ max.(a_eff_vec, 1e-12)
usd_pt          = usd_sh .> 0.0

# --- durations at every (t,i) ---
spells_periods = running_spells(adj)              # T×N
duration_years = vec(spells_periods) ./ per_year  # flatten

# ---- core moments over the entire panel ----
m1 = mean(_fin(duration_years))                       # mean duration
#m2 = mean(_fin(Float64.(usd_pt)))                     # USD participation rate

dwr_fin = _fin(d_wealth_ratio)
m3 = mean(dwr_fin)                                    # mean d/wealth
m4 = var(dwr_fin)                                     # var d/wealth
#m5 = var(_fin(Float64.(usd_pt)))                      # var USD participation

# ---- NEW: extra first-moment moments for identification ----

# Adjustment frequency: targets κ (and ties to sanity check)
adj_rate = mean(vec(adj) .> 0.0)
m6 = adj_rate

# Durable ownership share: targets ν and F^d
owner_share = mean(vec(d) .> 0.0)
m7 = owner_share

# Mean USD share among all (not just an indicator): targets F^t
#m8 = mean(_fin(usd_sh))

# Share of "heavily dollarized" agents: USD share > 0.5
#m9 = mean(vec(usd_sh) .> 0.5)

#dwr_usd = d_wealth_ratio[usd_pt .== 1.0]
#dwr_usd = _fin(dwr_usd)       # remove missing/NaN
#m10 = mean(dwr_usd)


#outmoms = [m1, m2, m3, m4, m5, m6, m7, m8, m9, m10]
outmoms = [m1,  m3, m4,  m6, m7]

if any(.!isfinite.(outmoms))
badix = findfirst(!isfinite, outmoms)
@warn "makemoments (panel) produced NaN/Inf" bad_moment=badix bad_value=outmoms[badix]
end

if settings.verbose
println("Moments over full panel (|r0|×N = $(length(r0))×$(size(adj,2)))")
println("m1 duration                     = ", m1)
println("m2 USD participation mean       = ", m2)
println("m3 d/wealth mean                = ", m3)
println("m4 d/wealth var                 = ", m4)
println("m5 USD participation var        = ", m5)
println("m6 adjustment rate              = ", m6)
println("m7 durable ownership share      = ", m7)
println("m8 mean USD share               = ", m8)
println("m9 share USD share > 0.5        = ", m9)
println("m10 mean d/wealth | USD holder  = ", m10)
end

return outmoms
end
