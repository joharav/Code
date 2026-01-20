@inline function is_col_stoch(P; tol=1e-10)
    maximum(abs.(sum(P, dims=1) .- 1.0)) < tol
end

@inline function is_row_stoch(P; tol=1e-10)
    maximum(abs.(sum(P, dims=2) .- 1.0)) < tol
end

"""
valfun(pea)

Wrapper that solves the adjust and no-adjust problems, then merges value + policies
using the value-maximizing regime choice.

Improvements vs your version:
- Enforces consistent array element types and shapes.
- Uses `@views` + preallocation to reduce allocations.
- Uses explicit `Bool` mask and avoids repeated `ifelse.` broadcasts where it matters.
- Optional transition-matrix sanity checks (row-stochastic for simulation draw code).
- Cleaner diagnostics and safer “share where adjust lower” computation.
"""
function valfun(pea::Vector{Float64}; grid_builder = makegrids)
    verbose = settings.verbose

    verbose && println("\n=== Solving ADJUSTMENT regime ===")
    adj = valfun_adjust(pea; grid_builder = makegrids)

    verbose && println("\n=== Solving NO-ADJUSTMENT regime ===")
    nadj = valfun_noadjust(pea; grid_builder = makegrids)

    # --- basic consistency checks (cheap, fail-fast) ---
    @assert size(adj.v) == size(nadj.v) "valfun: adj.v and nadj.v sizes differ"
    @assert size(adj.g.w) == size(nadj.g.w) "valfun: grids differ across regimes (should not happen)"

    # Optional: transition matrix sanity (your sim uses ROW-stochastic CDF rows)
    if hasproperty(adj.g, :t)
        P = adj.g.t
        if !(is_row_stoch(P) || is_col_stoch(P))
            @warn "Transition matrix is not row- or column-stochastic (check kron / tauchen)" maxerr_row=maximum(abs.(sum(P,dims=2).-1.0)) maxerr_col=maximum(abs.(sum(P,dims=1).-1.0))
        end
        if is_col_stoch(P) && !is_row_stoch(P)
            @warn "Transition matrix is column-stochastic but simulation expects row-stochastic; use P' in sim or build kron order accordingly."
        end
    end

    # --- regime choice mask ---
    # Use >= to break ties toward adjust (matches your original)
    @views choose_adj = adj.v .>= nadj.v              # Bool array
    adjust_flag = Float64.(choose_adj)                # keep if you use it in diagnostics

    # --- merged value ---
    v = similar(adj.v)
    @inbounds @simd for i in eachindex(v)
        v[i] = choose_adj[i] ? adj.v[i] : nadj.v[i]
    end

    # --- merged indices (Int arrays) ---
    # Allocate once, fill with loops (cheaper than multiple ifelse. broadcasts)
    iw = similar(adj.gidx.w)
    id = similar(adj.gidx.d)
    is = similar(adj.gidx.s)

    @inbounds for I in eachindex(iw)
        if choose_adj[I]
            iw[I] = adj.gidx.w[I]
            id[I] = adj.gidx.d[I]
            is[I] = adj.gidx.s[I]
        else
            iw[I] = nadj.gidx.w[I]
            id[I] = nadj.gidx.d[I]
            is[I] = nadj.gidx.s[I]
        end
    end
    gidx = dtp.Ipol(iw, id, is)

    # --- merged policies (Float64 arrays) ---
    pw  = similar(adj.pol.w)
    pd  = similar(adj.pol.d)
    ps  = similar(adj.pol.s)
    pc  = similar(adj.pol.c)
    paa = similar(adj.pol.aa)
    pa  = similar(adj.pol.a)

    @inbounds for I in eachindex(pw)
        if choose_adj[I]
            pw[I]  = adj.pol.w[I]
            pd[I]  = adj.pol.d[I]
            ps[I]  = adj.pol.s[I]
            pc[I]  = adj.pol.c[I]
            paa[I] = adj.pol.aa[I]
            pa[I]  = adj.pol.a[I]
        else
            pw[I]  = nadj.pol.w[I]
            pd[I]  = nadj.pol.d[I]
            ps[I]  = nadj.pol.s[I]
            pc[I]  = nadj.pol.c[I]
            paa[I] = nadj.pol.aa[I]
            pa[I]  = nadj.pol.a[I]
        end
    end
    pol = dtp.Pol(pw, pd, ps, pc, paa, pa)

    # --- diagnostics ---
    if verbose
        adj_prob = mean(adjust_flag)
        println("\nAdjustment probability (value-based): ", round(adj_prob, digits=3))

        # value difference stats without huge temporaries
        vdiff_mean = 0.0
        vdiff_min  = Inf
        vdiff_max  = -Inf
        nneg = 0
        n = length(adj.v)

        @inbounds for i in eachindex(adj.v)
            d = adj.v[i] - nadj.v[i]
            vdiff_mean += d
            vdiff_min = min(vdiff_min, d)
            vdiff_max = max(vdiff_max, d)
            nneg += (d < 0)
        end
        vdiff_mean /= n

        println("Mean diff: ", vdiff_mean, " | Min: ", vdiff_min, " | Max: ", vdiff_max)
        println("Share where adjust gives lower value: ", nneg / n)
        plotstuff4D(v,pol.w,pol.d,pol.s,pol.c,adj.g) 

   

    end

    # --- return bundle ---
    return (
        v = v,
        gidx = gidx,
        pol = pol,
        g = adj.g,                         # same grids in both
        e = max(adj.e, nadj.e),            # keep your stopping metric convention
        adjust_result = adj,
        noadjust_result = nadj,
        adjust_flag = adjust_flag,
        adjustment_indicator = choose_adj,
        pea = pea
    )




end
