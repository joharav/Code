# 4D Model Migration Guide

## Overview

This directory contains the refactored model code that reduces the state space from 5D to 4D:

**Original 5D Model:**
- States: `(e, y, aa, a, d)` - exchange rate, income, peso assets, dollar assets, durables
- Grid points: `11^5 = 161,051`
- Policy dimensions: `19 × 19 × 19 = 6,859` per state
- Total operations: ~1.96 billion per VFI iteration

**New 4D Model:**
- States: `(e, y, w, d)` - exchange rate, income, total wealth, durables
- Grid points: `7 × 7 × 15 × 11 = 8,085` (with finer grids)
- Policy dimensions: `21 × 21 × 15 = 6,615` per state (includes portfolio choice)
- Total operations: ~53 million per VFI iteration
- **Speedup: ~37x**

## Key Changes

### 1. State Space Transformation

The key insight is that `w = aa + e*a` (total wealth in pesos) is a sufficient state variable. The dollar share `s = e*a/w` becomes a **within-period choice** rather than a state.

### 2. Portfolio Choice

For each savings decision `w'`, the household optimally chooses:
```
s* = argmax E[V(e', y', w_realized, d')] - κ * s * w'
```
where:
```
w_realized = (1-s)*w'*(1+r) + s*w'*(1+r*)*(E'/E) - κ*s*w'
```

### 3. File Changes

| Original File | New File | Changes |
|---------------|----------|---------|
| `durable_mod.jl` | `durable_mod.jl` | New grid dimensions, `dtp.Ipol` and `dtp.Pol` structs now have `w`, `d`, `s` instead of `aa`, `a`, `d` |
| `makegrids.jl` | `makegrids.jl` | Single wealth grid `w` replaces `aa` and `a`; adds dollar share grid `s` |
| `utility.jl` | `utility.jl` | 6D utility array `[ie,iy,iw,id,iwp,idp]` |
| `utility_noadjust.jl` | `utility_noadjust.jl` | 5D utility array (no durable choice) |
| `maxbellman.jl` | `maxbellman.jl` | Integrated portfolio choice in inner loop |
| `maxbellman_noadjust.jl` | `maxbellman_noadjust.jl` | Same |
| `valfun_adjust.jl` | `valfun_adjust.jl` | 4D arrays throughout |
| `valfun_noadjust.jl` | `valfun_noadjust.jl` | 4D arrays throughout |
| `interpol.jl` | `interpol.jl` | Bilinear over `(w, d)` instead of trilinear |
| `simmodel.jl` | `simmodel.jl` | Tracks wealth evolution with portfolio returns |
| `makemoments.jl` | `makemoments.jl` | Derives `aa`, `a` from `w`, `s` for moments |
| `gmmfunctions.jl` | `gmmfunctions.jl` | Updated parameter bounds |
| `gridsearch_durables.jl` | `gridsearch_durables.jl` | New bounds for `F_d` and `κ` |

## Parameter Bound Changes

Based on diagnostic analysis showing severe moment mismatch:

| Parameter | Old Bounds | New Bounds | Rationale |
|-----------|------------|------------|-----------|
| `F_d` (durable fixed cost) | [0.01, 1.0] | [0.5, 5.0] | Model predicted 46% adjustment rate vs 8.6% in data |
| `κ` (dollar trans. cost) | [0.0, 1.0] | [0.01, 0.5] | Model predicted 12% dollar share vs 56% in data |

## Running the Model

### Test run:
```julia
include("collectfunctions.jl")
moms = test_4D_model()
```

### Estimation:
```julia
include("gridsearch_durables.jl")
```

## Files in this Directory

### Core Model (UPDATED FOR 4D)
- `durable_mod.jl` - Module definitions, grid sizes, data types
- `makegrids.jl` - Grid construction (w instead of aa, a)
- `tauchen.jl` - Tauchen discretization (unchanged)
- `interpol.jl` - Interpolation functions (bilinear over w, d)
- `ptrue.jl` - Default parameter values

### Utility (UPDATED)
- `utility.jl` - Flow utility (adjustment regime)
- `utility_noadjust.jl` - Flow utility (no adjustment)

### Value Function Iteration (UPDATED)
- `maxbellman.jl` - Full grid Bellman with portfolio choice (adjustment)
- `maxbellman_noadjust.jl` - Full grid Bellman (no adjustment)
- `tinybellman.jl` - Local search Bellman (adjustment)
- `tinybellman_noadjust.jl` - Local search Bellman (no adjustment)
- `howard.jl` - Howard acceleration (adjustment)
- `howard_noadjust.jl` - Howard acceleration (no adjustment)
- `valfun_adjust.jl` - VFI loop (adjustment)
- `valfun_noadjust.jl` - VFI loop (no adjustment)
- `valfun.jl` - Combined wrapper
- `fillin.jl` - Grid interpolation (state → policy)

### Policy & Simulation (UPDATED)
- `makepol.jl` - Policy extraction
- `simmodel.jl` - Panel simulation
- `simmodel_girf.jl` - Impulse response simulation
- `ergodic.jl` - Ergodic distribution

### Moments & Estimation (UPDATED)
- `makemoments.jl` - Moment computation
- `momentgen.jl` - Main moment generation wrapper
- `gmmfunctions.jl` - SMM estimation functions
- `gmmfunctions_broad.jl` - Broader SMM functions
- `gridsearch_durables.jl` - Estimation runner
- `simann.jl` - Simulated annealing (unchanged)
- `collectfunctions.jl` - Main include file

### Analysis (UPDATED)
- `diagnostics.jl` - Model diagnostics
- `diagnostics_runner.jl` - Diagnostics runner
- `counterfactuals.jl` - Counterfactual analysis
- `compstat.jl` - Comparative statics
- `decision_rules.jl` - Policy visualization

### Welfare (UPDATED)
- `welfare.jl` - Core welfare computations
- `welfare_cases.jl` - Welfare scenarios
- `welfare_disaster.jl` - Disaster welfare
- `welfare_dispersion.jl` - Welfare dispersion

### Disaster Analysis (UNCHANGED - may need updates)
- `disaster_grids.jl` - Disaster scenario grids
- `disaster_process.jl` - Disaster process
- `run_disaster_counterfactual.jl` - Disaster counterfactuals

### Plotting & Utilities (UNCHANGED)
- `plotgaps.jl`, `plotgaps_comp.jl`, `plotgaps_shock.jl`
- `plotstuff.jl`, `plotdensities.jl`
- `plot_distribution_panels.jl`, `plot_ergodic.jl`
- `heatmap.jl`
- `adj_gaps_sim.jl` - Adjustment gaps analysis
- `aggregate_series.jl`, `d_adjust_time_size.jl`
- `empirical_dist.jl`, `evalfun.jl`
- `girf.jl`, `inbetween.jl`, `mew.jl`
- `inflnc.jl`, `inflnc_functions.jl`
- `printstuff.jl`

## Computational Notes

1. **Memory:** The 4D model uses significantly less memory
2. **Parallelization:** All Bellman loops use `Threads.@threads`
3. **Interpolation:** Bilinear interpolation is faster than trilinear
4. **Portfolio choice:** Solved via grid search over `s` (15 points)

## Verification Checklist

- [ ] Test utility function produces reasonable consumption
- [ ] VFI converges for adjust and non-adjust regimes
- [ ] Simulation produces stationary distributions
- [ ] Moments are in reasonable ranges:
  - duration_mean: ~15-30 years
  - adj_rate: ~5-15%
  - dollar_share: ~40-70%
- [ ] Compare policy functions qualitatively to 5D model
