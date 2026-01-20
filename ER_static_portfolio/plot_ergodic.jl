using Plots

function plot_ergodic_w(dist4::Array{Float64,4}, g::NamedTuple, savepath::String)
    isdir(savepath) || mkpath(savepath)
    w = g.w
    ne, ny, nw, nd = size(dist4)

    for ie in 1:ne
        p = plot()
        for iy in 1:ny
            marg_w = sum(dist4[ie, iy, :, :], dims=2)[:]   # sum over d -> (nw,)
            s = sum(marg_w)
            if s > 0
                marg_w ./= s
                plot!(p, w, marg_w, label="iy=$iy")
            end
        end
        xlabel!(p, "Total wealth w")
        ylabel!(p, "Density")
        title!(p, "Ergodic wealth dist | e=$ie")
        savefig(p, joinpath(savepath, "w_dist_e$(ie).png"))
    end
    return nothing
end
function plot_ergodic_d_4D(dist4::Array{Float64,4}, g::NamedTuple, savepath::String)
    isdir(savepath) || mkpath(savepath)
    d = g.d
    ne, ny, nw, nd = size(dist4)

    for ie in 1:ne
        p = plot()
        for iy in 1:ny
            marg_d = sum(dist4[ie, iy, :, :], dims=1)[:]   # sum over w -> (nd,)
            s = sum(marg_d)
            if s > 0
                marg_d ./= s
                plot!(p, d, marg_d, label="iy=$iy")
            end
        end
        xlabel!(p, "Durables d")
        ylabel!(p, "Density")
        title!(p, "Ergodic durable dist | e=$ie")
        savefig(p, joinpath(savepath, "d_dist_e$(ie).png"))
    end
    return nothing
end
function plot_joint_heatmap_wd(dist4::Array{Float64,4}, g::NamedTuple, savepath::String)
    isdir(savepath) || mkpath(savepath)
    w, d = g.w, g.d
    ne, ny, nw, nd = size(dist4)

    for ie in 1:ne, iy in 1:ny
        H = dist4[ie, iy, :, :]             # (nw, nd)
        s = sum(H)
        if s > 0
            Hn = H ./ s
            # heatmap expects x, y and matrix with size (length(y), length(x)) if transposed
            p = heatmap(w, d, Hn', xlabel="Wealth w", ylabel="Durables d",
                        title="Joint (w,d) | e=$ie, y=$iy")
            savefig(p, joinpath(savepath, "joint_w_d_e$(ie)_y$(iy).png"))
        end
    end
    return nothing
end
function plot_total_marginals_4D(dist4::Array{Float64,4}, g::NamedTuple, savepath::String)
    isdir(savepath) || mkpath(savepath)
    w, d = g.w, g.d

    marg_w = sum(dist4, dims=(1,2,4))[:]   # (nw,)
    marg_d = sum(dist4, dims=(1,2,3))[:]   # (nd,)

    sw = sum(marg_w); sd = sum(marg_d)
    sw > 0 && (marg_w ./= sw)
    sd > 0 && (marg_d ./= sd)

    pW = plot(w, marg_w, xlabel="w", ylabel="Density", title="Aggregate wealth dist", legend=false)
    savefig(pW, joinpath(savepath, "total_w_dist.png"))

    pD = plot(d, marg_d, xlabel="d", ylabel="Density", title="Aggregate durable dist", legend=false)
    savefig(pD, joinpath(savepath, "total_d_dist.png"))

    return nothing
end
function plot_adjustment_by_income_4D(adjustment_indicator::AbstractArray{Bool,4},
    dist4::Array{Float64,4},
    savepath::String)
isdir(savepath) || mkpath(savepath)
ne, ny, nw, nd = size(dist4)

adj_by_income = fill(NaN, ny)
for iy in 1:ny
den = sum(dist4[:, iy, :, :])
if den > 0
num = sum(dist4[:, iy, :, :][adjustment_indicator[:, iy, :, :]])
adj_by_income[iy] = num / den
end
end

p = plot(1:ny, adj_by_income, marker=:circle,
xlabel="Income state", ylabel="Adjustment rate",
title="Adjustment rate by income", legend=false)
savefig(p, joinpath(savepath, "adjustment_by_income.png"))
return nothing
end
