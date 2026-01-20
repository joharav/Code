using StatsBase

# --- helpers ---
@inline function idx_on_grid(x::Real, grid::AbstractVector{<:Real})
    # returns nearest lower bracket index in 1..length(grid)
    j = searchsortedlast(grid, x)
    return clamp(j, 1, length(grid))
end

@inline function clamp_level(x::Real, grid::AbstractVector{<:Real})
    return clamp(x, grid[1], grid[end])
end

"""
    empirical_init_4D(nagents, df, grids;
                      w_col=:w, d_col=:d, y_col=:y, ex_col=:ex,
                      a_usd_col=:a, aa_loc_col=:aa, peso_col=:pesoEFHU,
                      build_w_from_assets=true, default_ex=:mean)

Return initial indices for current setting (e,y,w,d) and (optionally) s.

- If `build_w_from_assets=true`, we compute w and s from aa_loc + ex*a_usd.
- If you already have total liquid wealth `w_col` and dollar share `s_col`,
  set `build_w_from_assets=false` and provide those columns.
"""
function empirical_init_4D(nagents::Int,
                           df::DataFrame,
                           grids::NamedTuple;
                           # column names (set these to your EFHU merged names)
                           d_col::Symbol = :durables_final,
                           y_col::Symbol = :household_income_q,
                           ex_col::Union{Symbol,Nothing} = :exrate,   # set to nothing if not in df
                           a_usd_col::Symbol = :a_usd,                # USD assets (level)
                           aa_loc_col::Symbol = :aa_loc,              # local assets (level)
                           w_col::Symbol = :w,                        # only used if build_w_from_assets=false
                           s_col::Symbol = :s,                        # only used if build_w_from_assets=false
                           peso_col::Union{Symbol,Nothing} = :pesoEFHU,
                           build_w_from_assets::Bool = true,
                           default_ex::Symbol = :mean)

    # --- weights ---
    if peso_col === nothing || !(peso_col in names(df))
        weights = Weights(fill(1.0, nrow(df)))
    else
        wts = df[!, peso_col]
        wts = float.(coalesce.(wts, 0.0))
        # guard: if weights are degenerate, fall back to uniform
        if sum(wts) <= 0
            weights = Weights(fill(1.0, nrow(df)))
        else
            weights = Weights(wts ./ sum(wts))
        end
    end

    # --- draw micro sample ---
    draw = sample(1:nrow(df), nagents; replace=true, weights=weights)
    sdf = df[draw, :]

    # --- exchange rate levels for initialization ---
    ex_level = Vector{Float64}(undef, nagents)
    if ex_col !== nothing && (ex_col in names(sdf))
        ex_level .= float.(sdf[!, ex_col])
        # clamp to grid support
        @inbounds for i in 1:nagents
            ex_level[i] = clamp_level(ex_level[i], grids.ex)
        end
    else
        # no exchange-rate column in data: initialize at some default point
        ex0 = default_ex === :mean ? mean(grids.ex) :
              default_ex === :mid  ? grids.ex[cld(length(grids.ex),2)] :
              grids.ex[1]
        ex0 = clamp_level(ex0, grids.ex)
        fill!(ex_level, ex0)
    end

    # --- build d and y levels ---
    d_level = clamp_level.(float.(sdf[!, d_col]), Ref(grids.d))
    y_level = clamp_level.(float.(sdf[!, y_col]), Ref(grids.y))

    # --- build w and s levels (current setting) ---
    w_level = Vector{Float64}(undef, nagents)
    s_level = Vector{Float64}(undef, nagents)

    if build_w_from_assets
        a_usd = float.(sdf[!, a_usd_col])
        aa_loc = float.(sdf[!, aa_loc_col])

        @inbounds for i in 1:nagents
            a_i  = max(a_usd[i], 0.0)
            aa_i = max(aa_loc[i], 0.0)
            e_i  = ex_level[i]

            w_i = aa_i + e_i * a_i
            w_i = clamp_level(w_i, grids.w)

            s_i = w_i > 0 ? (e_i * a_i) / w_i : 0.0
            s_i = clamp(s_i, 0.0, 1.0)

            w_level[i] = w_i
            s_level[i] = s_i
        end
    else
        w_level .= clamp_level.(float.(sdf[!, w_col]), Ref(grids.w))
        s_level .= clamp.(float.(sdf[!, s_col]), 0.0, 1.0)
    end

    # --- map to grid INDICES (this is what your simulator should store) ---
    ie0 = [idx_on_grid(ex_level[i], grids.ex) for i in 1:nagents]
    iy0 = [idx_on_grid(y_level[i],  grids.y)  for i in 1:nagents]
    iw0 = [idx_on_grid(w_level[i],  grids.w)  for i in 1:nagents]
    id0 = [idx_on_grid(d_level[i],  grids.d)  for i in 1:nagents]
    is0 = [idx_on_grid(s_level[i],  grids.s)  for i in 1:nagents]  # if you have an s grid

    return (
        # indices
        ie = ie0,
        iy = iy0,
        iw = iw0,
        id = id0,
        is = is0,
        # optional levels for debugging
        ex = ex_level,
        y  = y_level,
        w  = w_level,
        d  = d_level,
        s  = s_level
    )
end
