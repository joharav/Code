using LinearAlgebra, Statistics

# =======================
# Basic influence functions
# =======================

"IF of mean(x)"
mean_if(x::AbstractVector{<:Real}) = x .- mean(x)

"IF of variance(x) around the sample mean"
function var_if(x::AbstractVector{<:Real})
    μ = mean(x)
    v = mean((x .- μ).^2)
    return (x .- μ).^2 .- v
end

"IF of covariance(x,y) around sample means"
function cov_if(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y)
    mx, my = mean(x), mean(y)
    c = mean((x .- mx) .* (y .- my))
    return (x .- mx) .* (y .- my) .- c
end

"IF of a share/probability: θ = mean(1{x>0})"
function share_if(x::AbstractVector{<:Real})
    θ = mean(x .> 0.0)
    return (x .> 0.0) .- θ
end

"IF of a ratio of means: θ = mean(x) / mean(y)"
function ratio_of_means_if(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    @assert length(x) == length(y)
    mx = mean(x); my = mean(y)
    return (x .- mx)./my .- (mx/(my^2)) .* (y .- my)
end

# =======================
# OLS (coefficients) + IF
# =======================

"OLS coefficients via a stable solver"
ols_coef(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}) = X \ y

"""
IF for OLS β̂ (rows = observations, columns = parameters).

Formula: IF_i(β) = (X'X / n)^(-1) * x_i * u_i
"""
function ols_if(y::AbstractVector{<:Real}, X::AbstractMatrix{<:Real})
    n, k = size(X)
    @assert length(y) == n
    β = ols_coef(X, y)
    u = y .- X*β
    # (X'X/n)^(-1) = n * (X'X)^{-1}; we build it via a factorization for stability
    XtX = Symmetric(X'X)
    F = cholesky(XtX)             # SPD assumed; use qr if X not full rank
    XtX_inv = inv(F) * inv(F)'    # same as (X'X)^{-1}
    A = n * XtX_inv               # (X'X/n)^(-1)
    # IF_i = A * (x_i * u_i)  -> stack all obs into a matrix of size n×k
    return (X .* u) * A'          # n×k
end

# =======================
# (Double-)Clustered covariance from IFs
# =======================

"""
CRV1-style cluster covariance from influence functions.

Inputs:
  φ : n×p matrix of influence functions (row i = IF_i(θ)', p params)
  g : length-n vector of cluster ids

Returns: p×p covariance matrix ≈ (∑_c s_c s_c') / n^2
where s_c = ∑_{i∈c} φ_i
"""
function cluster_cov(φ::AbstractMatrix{<:Real}, g)
    n, p = size(φ)
    gids = unique(g)
    V = zeros(p, p)
    for id in gids
        phic = φ[g .== id, :]
        sc = vec(sum(phic, dims=1))
        V .+= sc * sc'
    end
    return V / n^2
end

"""
Cameron–Gelbach–Miller double clustering:

V = V_g1 + V_g2 - V_naive

where V_naive = (φ'φ)/n^2.
"""
function double_cluster_cov(φ::AbstractMatrix{<:Real}, g1, g2)
    n = size(φ,1)
    V1 = cluster_cov(φ, g1)
    V2 = cluster_cov(φ, g2)
    Vn = (φ'φ) / n^2
    return V1 .+ V2 .- Vn
end

# =======================
# Weighted helpers (optional)
# =======================

"Weighted mean that tolerates missings"
function wmean(x::AbstractVector{<:Union{Missing,Real}}, w::AbstractVector{<:Real})
    v = .!ismissing.(x)
    xv = Float64.(x[v]); wv = w[v]
    return sum(wv .* xv) / sum(wv)
end

"Weighted variance around weighted mean (no small-sample correction)"
function wvar(x::AbstractVector{<:Union{Missing,Real}}, w::AbstractVector{<:Real})
    v = .!ismissing.(x)
    xv = Float64.(x[v]); wv = w[v]
    μ = sum(wv .* xv) / sum(wv)
    return sum(wv .* (xv .- μ).^2) / sum(wv)
end
