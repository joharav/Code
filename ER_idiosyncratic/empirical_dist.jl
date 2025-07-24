using CSV, DataFrames

df = CSV.read("EFHU_init.csv", DataFrame)

# Optional: sample according to pesoEFHU
weights = df.pesoEFHU ./ sum(df.pesoEFHU)
sample_ids = sample(1:nrow(df), sz.nagents; replace=true, weights=weights)

sampled = df[sample_ids, :]

# Discretize to match grids
function clamp_to_grid(val, grid)
    return grid[clamp(searchsortedfirst(grid, val), 1, length(grid))]
end

a_empirical = [clamp_to_grid(a, g.a) for a in sampled.total_savings]
d_empirical = [clamp_to_grid(d, g.d) for d in sampled.durables]
y_empirical = [clamp_to_grid(y, g.y) for y in sampled.household_income]
ex_empirical = fill(mean(g.ex), sz.nagents) # Or a default value

# Return as tuple or NamedTuple
init_state = (
    a = a_empirical,
    d = d_empirical,
    y = y_empirical,
    ex = ex_empirical
)
function empirical_init(nagents::Int, df::DataFrame, grids::NamedTuple)
    # Extract relevant variables from EFHU data frame
    d_raw = df.d
    a_raw = df.a
    y_raw = df.y
    e_raw = df.e  # if available, or generate synthetic exchange shocks

    # Optional: apply transformation or clamp to model grid support
    d_grid_idx = [searchsortedfirst(grids.d, d) for d in d_raw]
    a_grid_idx = [searchsortedfirst(grids.a, a) for a in a_raw]
    y_idx = [searchsortedfirst(grids.y, y) for y in y_raw]
    e_idx = [searchsortedfirst(grids.ex, e) for e in e_raw]

    return (
        a = a_grid_idx,
        d = d_grid_idx,
        y = y_idx,
        ex = e_idx
    )
end
