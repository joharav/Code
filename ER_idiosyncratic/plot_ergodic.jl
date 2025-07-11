using Plots

function plot_ergodic_a(dist::Array{Float64,4}, g::NamedTuple)
    a = g.a
    ne, ny, nd, na = size(dist)

    for ie in 1:ne
        plot()
        for iy in 1:ny
            # Sum over durables, normalize
            marg_a = sum(dist[ie, iy, :, :], dims=2)[:]
            marg_a ./= sum(marg_a)

            plot!(a, marg_a, label="iy=$iy")
        end
        xlabel!("Assets")
        ylabel!("Density")
        title!("Ergodic Asset Distribution | Exchange rate state $ie")
        savefig("Output/Ergodic/a_dist_e$(ie).png")
    end
end
function plot_ergodic_d(dist::Array{Float64,4}, g::NamedTuple)
    d = g.d
    ne, ny, nd, na = size(dist)

    for ie in 1:ne
        plot()
        for iy in 1:ny
            # Sum over assets, normalize
            marg_d = sum(dist[ie, iy, :, :], dims=1)[:]
            marg_d ./= sum(marg_d)

            plot!(d, marg_d, label="iy=$iy")
        end
        xlabel!("Durables")
        ylabel!("Density")
        title!("Ergodic Durable Distribution | Exchange rate state $ie")
        savefig("Output/Ergodic/d_dist_e$(ie).png")
    end
end
function plot_joint_heatmap(dist::Array{Float64,4}, g::NamedTuple)
    d = g.d
    a = g.a
    ne, ny = sz.ne, sz.ny

    for ie in 1:ne, iy in 1:ny
        heat_data = dist[ie, iy, :, :] ./ sum(dist[ie, iy, :, :])
        heatmap(a, d, heat_data', xlabel="Assets", ylabel="Durables",
                title="Joint (a,d) | e=$ie, y=$iy", colorbar_title="Density")
        savefig("Output/Ergodic/joint_a_d_e$(ie)_y$(iy).png")
    end
end
