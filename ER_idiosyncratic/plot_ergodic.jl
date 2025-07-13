using Plots

function plot_ergodic_a(dist::Array{Float64,4}, g::NamedTuple, savepath::String)
    a = g.a
    ne, ny, nd, na = size(dist)

    for ie in 1:ne
        plot()
        for iy in 1:ny
            marg_a = sum(dist[ie, iy, :, :], dims=2)[:]
            marg_a ./= sum(marg_a)
            plot!(a, marg_a, label="iy=$iy")
        end
        xlabel!("Assets")
        ylabel!("Density")
        title!("Ergodic Asset Distribution | Exchange rate state $ie")
        savefig("$savepath/a_dist_e$(ie).png")
    end
end
function plot_ergodic_d(dist::Array{Float64,4}, g::NamedTuple, savepath::String)
    d = g.d
    ne, ny, nd, na = size(dist)

    for ie in 1:ne
        plot()
        for iy in 1:ny
            marg_d = sum(dist[ie, iy, :, :], dims=1)[:]
            marg_d ./= sum(marg_d)
            plot!(d, marg_d, label="iy=$iy")
        end
        xlabel!("Durables")
        ylabel!("Density")
        title!("Ergodic Durable Distribution | Exchange rate state $ie")
        savefig("$savepath/d_dist_e$(ie).png")
    end
end

function plot_joint_heatmap(dist::Array{Float64,4}, g::NamedTuple, savepath::String)
    d = g.d
    a = g.a
    ne, ny = sz.ne, sz.ny

    for ie in 1:ne, iy in 1:ny
        heat_data = dist[ie, iy, :, :] ./ sum(dist[ie, iy, :, :])
        heatmap(a, d, heat_data', xlabel="Assets", ylabel="Durables",
                title="Joint (a,d) | e=$ie, y=$iy", colorbar_title="Density")
        savefig("$savepath/joint_a_d_e$(ie)_y$(iy).png")
    end
end

function plot_adjustment_by_income(adjustment_indicator::BitArray{4}, dist::Array{Float64,4}, savepath::String)
    ne, ny, nd, na = size(dist)
    adj_by_income = zeros(ny)

    for iy in 1:ny
        numerator = sum(dist[:, iy, :, :][adjustment_indicator[:, iy, :, :]])
        denominator = sum(dist[:, iy, :, :])
        adj_by_income[iy] = numerator / denominator
    end

    plot(1:ny, adj_by_income, marker=:circle, xlabel="Income state", ylabel="Adjustment Rate",
        title="Adjustment Rate by Income Level", legend=false)
    savefig("$savepath/adjustment_by_income.png")
end



function plot_total_marginals(dist::Array{Float64,4}, g::NamedTuple, savepath::String)
    a, d = g.a, g.d

    marg_a = sum(dist, dims=(1, 2, 3))[:]
    marg_a ./= sum(marg_a)
    plot(a, marg_a, xlabel="Assets", ylabel="Density", title="Aggregate Asset Distribution", legend=false)
    savefig("$savepath/total_asset_dist.png")

    marg_d = sum(dist, dims=(1, 2, 4))[:]
    marg_d ./= sum(marg_d)
    plot(d, marg_d, xlabel="Durables", ylabel="Density", title="Aggregate Durable Distribution", legend=false)
    savefig("$savepath/total_durable_dist.png")
end


function plot_cev_distribution(cev::Array{Float64,4}, dist::Array{Float64,4})
    cev_vec = cev[:]
    weights = dist[:]
    weights ./= sum(weights)

    # Weighted histogram
    histogram(cev_vec, weights=weights, bins=50, xlabel="CEV", ylabel="Density",
        title="Welfare Gains (CEV) from Currency Regime", legend=false)
    savefig("Output/Ergodic/cev_distribution.png")
end