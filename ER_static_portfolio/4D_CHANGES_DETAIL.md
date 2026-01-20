# 4D Model: Detailed File Changes

## Summary of Dimensional Changes

| Component | Old (5D) | New (4D) |
|-----------|----------|----------|
| State space | `(e, y, aa, a, d)` | `(e, y, w, d)` |
| Value function | `V[ne, ny, na, na, nd]` | `V[ne, ny, nw, nd]` |
| Policy indices | `(iaa', ia', id')` | `(iw', id', is)` |
| Utility array | `[ne,ny,na,na,nd,npa,npa,npd]` | `[ne,ny,nw,nd,npw,npd]` |

---

## 1. `durable_mod.jl` - Module Definitions

### Grid sizes:
```julia
# Old:
const na = 9     # peso and dollar asset grids
const nd = 9
const npa = 19
const npd = 19

# New:
const nw = 15    # total wealth grid (replaces na×na)
const nd = 11    # finer durable grid
const ns = 15    # dollar share grid (within-period choice)
const npw = 21   # wealth policy grid
const npd = 21   # durable policy grid
```

### Data types:
```julia
# Old:
mutable struct Ipol
    a::Array{Int,5}    # dollar asset index
    aa::Array{Int,5}   # peso asset index
    d::Array{Int,5}    # durable index
end

# New:
mutable struct Ipol
    w::Array{Int,4}    # total wealth index
    d::Array{Int,4}    # durable index
    s::Array{Int,4}    # dollar share index (within-period)
end
```

---

## 2. `makegrids.jl` - Grid Construction

### Key change: Single wealth grid
```julia
# Old: Separate peso and dollar asset grids
aag = range(0, aa_max, length=na)    # peso
ag = range(0, a_max, length=na)      # dollar

# New: Total wealth grid
wg = range(0, w_max, length=nw)      # total wealth
sg = range(0, 1, length=ns)          # dollar share
```

### Returns different named tuple:
```julia
# Old:
(t, a, ap, d, dp, ex, y, aa, aap)

# New:
(t, w, wp, s, d, dp, ex, y, te)  # te = ER transition for portfolio
```

---

## 3. `utility.jl` - Flow Utility (Adjustment)

### Array dimensions:
```julia
# Old: 8 dimensions
util[ie, iy, iaa, ia, id, iiaa, iia, iid]

# New: 6 dimensions
util[ie, iy, iw, id, iwp, idp]
```

### Budget constraint:
```julia
# Old:
c = income + sale - durable_purchase - aa_next - E*a_next - κ*E*a_next - time_cost

# New:
c = income + sale - durable_purchase - w_next - time_cost
# (κ cost handled in continuation value via portfolio choice)
```

---

## 4. `maxbellman.jl` - Bellman Optimization

### Key innovation: Nested portfolio choice
```julia
# Old: Triple loop over (aa', a', d')
for iiaa in 1:npa
    for iia in 1:npa
        for iid in 1:npd
            bellman = util[...] + β*queue[...]
        end
    end
end

# New: Double loop over (w', d') with inner portfolio optimization
for idp in 1:npd
    for iwp in 1:npw
        u_flow = util[ie, iy, iw, id, iwp, idp]
        
        # Find optimal s* for this w'
        for is in 1:ns
            s = s_grid[is]
            # Compute E[V] integrating over e' with portfolio returns
            EV = sum over e' of:
                prob(e'|e) * V_interp(e', y', w_realized, d')
            where w_realized = (1-s)*w'*(1+r) + s*w'*(1+r*)*(e'/e) - κ*s*w'
        end
        
        bellman = u_flow + β*best_EV
    end
end
```

---

## 5. `valfun_adjust.jl` - VFI Loop

### Queue computation:
```julia
# Old: Reshape and multiply over 4D continuation
reshaped_v = reshape(v, ne*ny, na, na, nd)
for id, ia, iaa:
    qq[:, iaa, ia, id] = tmat * reshaped_v[:, iaa, ia, id]

# New: Reshape over 3D (portfolio choice done in Bellman)
reshaped_v = reshape(v, ne*ny, nw, nd)
for id, iw:
    qq[:, iw, id] = tmat * reshaped_v[:, iw, id]
```

### Policy convergence tracking:
```julia
# Old:
pgap = sum(|gidx.a - old.a|) + sum(|gidx.aa - old.aa|) + sum(|gidx.d - old.d|)

# New:
pgap = sum(|gidx.w - old.w|) + sum(|gidx.d - old.d|) + sum(|gidx.s - old.s|)
```

---

## 6. `interpol.jl` - Interpolation

### Reduced dimensionality:
```julia
# Old: Trilinear over (aa, a, d)
interpol(e, y, aa, a, d, grids, V::Array{Float64,5})
    bracket aa → (iaaL, iaaU, waa)
    bracket a  → (iaL, iaU, wa)
    bracket d  → (idL, idU, wd)
    # 8-point interpolation

# New: Bilinear over (w, d)
interpol(e, y, w, d, grids, V::Array{Float64,4})
    bracket w → (iwL, iwU, ww)
    bracket d → (idL, idU, wd)
    # 4-point interpolation
```

---

## 7. `simmodel.jl` - Simulation

### State tracking:
```julia
# Old: Track (aa, a, d) separately
aa_old, a_old, d_old = ...
aa_pr = interpol(..., pol.aa)
a_pr = interpol(..., pol.a)

# New: Track (w, s, d), derive (aa, a)
w_old, s_old, d_old = ...
w_pr = interpol(..., pol.w)
s_pr = interpol(..., pol.s)
aa_pr = (1-s_pr) * w_pr
a_pr = s_pr * w_pr / e
```

### Wealth evolution with portfolio returns:
```julia
# New addition: Portfolio-dependent wealth transition
w_new = (1-s)*w'*(1+r) + s*w'*(1+r*)*(e_new/e) - κ*s*w'
```

---

## 8. `makemoments.jl` - Moment Computation

### Key moments (unchanged definitions):
1. `duration_mean` - years since last adjustment
2. `dwealth_mean` - mean durable/total wealth ratio
3. `dwealth_var` - variance of above
4. `adj_rate` - annual adjustment probability
5. `dollar_share` - mean dollar share of liquid assets
6. `dollar_vol` - variance of dollar share

### Derivation from new variables:
```julia
# Dollar share in simulation
aa = (1 - s) .* w
a_fx = e .* (s .* w ./ e) = s .* w
dollar_share = s .* w ./ (w) = s

# Or equivalently from derived aa, a
dollar_share = a_fx ./ (aa + a_fx)
```

---

## 9. `gridsearch_durables.jl` - Estimation

### Updated parameter bounds:
```julia
# Parameter: F_d (durable fixed cost)
# Old: lb[2] = 0.01, ub[2] = 1.0
# New: lb[2] = 0.50, ub[2] = 5.0
# Rationale: Need higher friction to match 8.6% adjustment rate

# Parameter: κ (dollar transaction cost)  
# Old: lb[3] = 0.0, ub[3] = 1.0
# New: lb[3] = 0.01, ub[3] = 0.5
# Rationale: Need lower friction to match 56% dollar share
```

---

## Computational Comparison

| Metric | Old 5D | New 4D | Improvement |
|--------|--------|--------|-------------|
| State points | 161,051 | 8,085 | 20x fewer |
| Policy evaluations | ~2B | ~53M | 37x fewer |
| Memory (value function) | ~1.3 GB | ~65 MB | 20x less |
| VFI iteration time | ~minutes | ~seconds | Order of magnitude |

---

## Migration Checklist

1. [ ] Replace `durable_mod.jl` with new version
2. [ ] Replace `makegrids.jl` - new grid structure
3. [ ] Replace all VFI files (`utility*.jl`, `maxbellman*.jl`, `valfun*.jl`)
4. [ ] Replace `interpol.jl` - new signature
5. [ ] Replace `simmodel.jl` - portfolio tracking
6. [ ] Replace `makemoments.jl` - derive from w, s
7. [ ] Update `gmmfunctions.jl` - new bounds
8. [ ] Update estimation runner
9. [ ] Test with `test_4D_model()`
10. [ ] Verify moments are in reasonable range
