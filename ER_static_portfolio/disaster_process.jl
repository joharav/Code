#-------------------------------------------------------------------------------
# Refactored Julia Code
#-------------------------------------------------------------------------------

using Random, Distributions, LinearAlgebra

# --- Core Functions ----------------------------------------------------------

"""
    annual_to_period_prob(pi_y::Float64, periods_per_year::Int)

Maps an annual probability `pi_y` to its equivalent per-period probability,
assuming independent trials within the year.

The formula is `1 - (1 - pi_y)^(1 / periods_per_year)`.

# Arguments
- `pi_y::Float64`: The annual probability, must be in `[0, 1]`.
- `periods_per_year::Int`: The number of periods in a year (e.g., 4 for quarterly).

# Returns
- `Float64`: The per-period probability.
"""
function annual_to_period_prob(pi_y::Float64, periods_per_year::Int)
    # This function is already great. The main change is adopting the standard
    # Julia docstring format for better integration with help systems (e.g., ? operator).
    @assert 0.0 ≤ pi_y ≤ 1.0 "pi_y must be in [0,1]"
    @assert periods_per_year ≥ 1 "periods_per_year must be ≥ 1"
    return 1.0 - (1.0 - pi_y)^(1 / periods_per_year)
end

"""
    EV_mixture_over_jump(loge, params, eps_nodes, eps_weights, next_state_iter)

Computes the conditional expectation `E[V(s') | s]` over a rare jump in the
exchange rate process.

This calculates `(1-πe) * E[V | no jump] + πe * E[V | jump]`, where the
expectation over the continuous shock `ε` is handled by quadrature.

# Arguments
- `loge::Float64`: Current log exchange rate.
- `params::NamedTuple`: A tuple containing process parameters (`ρ`, `σ`, `κe`, `πe`).
- `eps_nodes::AbstractVector`: Quadrature nodes for the `N(0,1)` innovation.
- `eps_weights::AbstractVector`: Quadrature weights for the `N(0,1)` innovation.
- `next_state_iter::Function`: A function `(loge_next) -> Float64` that computes
  the value function integrated over all *other* state variables for a given
  next-period log exchange rate `loge_next`.

# Returns
- `Float64`: The expected value.
"""
function EV_mixture_over_jump(
    loge::Float64,
    params::NamedTuple, # Suggestion: Group parameters for cleaner function signatures.
    eps_nodes::AbstractVector,
    eps_weights::AbstractVector,
    next_state_iter::F # Using a type parameter F is slightly more robust for closures.
) where {F<:Function}
    # Unpack parameters
    ρ, σ, κe, πe = params.ρ, params.σ, params.κe, params.πe

    EV_nojump = 0.0
    EV_jump = 0.0
    
    # The core logic is efficient and correct. No changes needed here.
    # Using @inbounds is fine as long as you guarantee nodes/weights match length.
    @inbounds for (ε, w) in zip(eps_nodes, eps_weights)
        loge_next_nojump = ρ * loge + σ * ε
        
        # Integrate over other states (e.g., income) via the provided closure
        val_nojump = next_state_iter(loge_next_nojump)
        val_jump = next_state_iter(loge_next_nojump + κe)

        EV_nojump += w * val_nojump
        EV_jump += w * val_jump
    end

    return (1.0 - πe) * EV_nojump + πe * EV_jump
end


# --- Example Usage ----------------------------------------------------------

# Suggestion: Group model parameters into a struct or NamedTuple. This makes
# passing them to functions much cleaner and less error-prone.
params_e = (
    ρ = 0.66,
    σ = 0.25,
    πe_annual = 0.02,
    periods_per_year = 4,
    κe = log(1.4)
)

# Derive the per-period probability from the annual one
πe = annual_to_period_prob(params_e.πe_annual, params_e.periods_per_year)

# Add it to the parameters tuple
params_e = merge(params_e, (πe = πe,))


# --- 1) Quadrature nodes/weights for ε ~ N(0,1) ---
# Your method using Gauss-Hermite is standard and correct.
function gauss_hermite_for_standard_normal(z, w)
    nodes = sqrt(2.0) .* z
    weights = (1.0 / sqrt(π)) .* w
    return nodes, weights
end

const GH_Z = [-2.020182870, -0.958572465, 0.0, 0.958572465, 2.020182870]
const GH_W = [0.0199532421, 0.393619323, 0.945308720, 0.393619323, 0.0199532421]
eps_nodes_e, eps_weights_e = gauss_hermite_for_standard_normal(GH_Z, GH_W)


# --- 2) Toy income Markov chain ---
y_trans_probs = [0.5, 0.5]
y_next_indices = [1, 2]


# --- 3) Fake V interpolator ---
function V_interp(a_p::Float64, aD_p::Float64, d_p::Float64, loge_p::Float64, yidx_p::Int)
    base = -(loge_p - 0.0)^2
    ybonus = (yidx_p == 2 ? 0.25 : 0.0)
    return base + ybonus
end

# --- 4) Choices and current state for the test ---
a_p, aD_p, d_p = 1.0, 0.5, 1.2
loge = log(1.0)

# Build the closure
# This is a great way to structure the problem.
next_state_iter = function (loge_next::Float64)
    EV_y = 0.0
    for (yidxp, Py) in zip(y_next_indices, y_trans_probs)
        EV_y += Py * V_interp(a_p, aD_p, d_p, loge_next, yidxp)
    end
    return EV_y
end

# --- 5) Compute the mixture expectation ---
EV = EV_mixture_over_jump(
        loge,
        params_e, # Pass the single parameters object
        eps_nodes_e,
        eps_weights_e,
        next_state_iter
)
println("EV (quadrature) = ", EV)
