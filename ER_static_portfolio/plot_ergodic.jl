using Plots

# --- asset (foreign 'a') marginal, by e-state ---
function plot_ergodic_a(dist5::Array{Float64,5}, g::NamedTuple, savepath::String)
    a = g.a
    ne, ny, nd, naa, na = size(dist5)
    for ie in 1:ne
        plot()
        for iy in 1:ny
            # sum over (id, iaa) to leave ia (foreign assets)
            marg_a = sum(dist5[ie, iy, :, :, :], dims=(3,4))[:]
            marg_a ./= max(sum(marg_a), eps())
            plot!(a, marg_a, label="iy=$iy")
        end
        xlabel!("Foreign assets a")
        ylabel!("Density")
        title!("Ergodic foreign-asset dist | e=$ie")
        savefig("$savepath/a_dist_e$(ie).png")
    end
end

# --- asset (local 'aa') marginal, by e-state ---
function plot_ergodic_aa(dist5::Array{Float64,5}, g::NamedTuple, savepath::String)
    aa = g.aa
    ne, ny, nd, naa, na = size(dist5)
    for ie in 1:ne
        plot()
        for iy in 1:ny
            # sum over (id, ia) to leave iaa (local assets)
            marg_aa = sum(dist5[ie, iy, :, :, :], dims=(3,5))
            v = vec(marg_aa)
            v ./= max(sum(v), eps())
            plot!(aa, v, label="iy=$iy")
        end
        xlabel!("Local assets aa")
        ylabel!("Density")
        title!("Ergodic local-asset dist | e=$ie")
        savefig("$savepath/aa_dist_e$(ie).png")
    end
end

# --- durables marginal, by e-state (unchanged logic, 5D input) ---
function plot_ergodic_d(dist5::Array{Float64,5}, g::NamedTuple, savepath::String)
    d = g.d
    ne, ny, nd, naa, na = size(dist5)
    for ie in 1:ne
        plot()
        for iy in 1:ny
            marg_d = sum(dist5[ie, iy, :, :, :], dims=(4,5))[:]
            marg_d ./= max(sum(marg_d), eps())
            plot!(d, marg_d, label="iy=$iy")
        end
        xlabel!("Durables d")
        ylabel!("Density")
        title!("Ergodic durable dist | e=$ie")
        savefig("$savepath/d_dist_e$(ie).png")
    end
end

# --- joint heatmap of (a,d) after summing out local aa ---
function plot_joint_heatmap(dist5::Array{Float64,5}, g::NamedTuple, savepath::String)
    d, a = g.d, g.a
    ne, ny = sz.ne, sz.ny
    for ie in 1:ne, iy in 1:ny
        # sum over iaa, keep (id, ia)
        M = sum(dist5[ie, iy, :, :, :], dims=4)  # => (id,1,ia)
        H = dropdims(M, dims=2)                  # (id, ia)
        H ./= max(sum(H), eps())
        heatmap(a, d, H', xlabel="Foreign assets a", ylabel="Durables d",
                title="Joint (a,d) | e=$ie, y=$iy", colorbar_title="Density")
        savefig("$savepath/joint_a_d_e$(ie)_y$(iy).png")
    end
end

# --- total marginals for (a, aa, d) ---
function plot_total_marginals(dist5::Array{Float64,5}, g::NamedTuple, savepath::String)
    a, aa, d = g.a, g.aa, g.d

    marg_a  = sum(dist5, dims=(1,2,3,4))[:] ; marg_a  ./= max(sum(marg_a), eps())
    marg_aa = sum(dist5, dims=(1,2,3,5))[:] ; marg_aa ./= max(sum(marg_aa), eps())
    marg_d  = sum(dist5, dims=(1,2,4,5))[:] ; marg_d  ./= max(sum(marg_d), eps())

    plot(a, marg_a,  xlabel="a",  ylabel="Density", title="Aggregate foreign-asset dist", legend=false)
    savefig("$savepath/total_a_dist.png")

    plot(aa, marg_aa, xlabel="aa", ylabel="Density", title="Aggregate local-asset dist", legend=false)
    savefig("$savepath/total_aa_dist.png")

    plot(d, marg_d,  xlabel="d",  ylabel="Density", title="Aggregate durable dist", legend=false)
    savefig("$savepath/total_d_dist.png")
end


function plot_adjustment_by_income(adjustment_indicator::BitArray{5},
    dist5::Array{Float64,5},
    savepath::String)
    _, ny, _, _, _ = size(dist5)
    adj_by_income = zeros(ny)
    for iy in 1:ny
        num = sum(dist5[:, iy, :, :, :][adjustment_indicator[:, iy, :, :, :]])
        den = sum(dist5[:, iy, :, :, :])
        adj_by_income[iy] = den > 0 ? num / den : NaN
    end
    plot(1:ny, adj_by_income, marker=:circle, xlabel="Income state",
    ylabel="Adjustment Rate", title="Adjustment Rate by Income Level", legend=false)
    savefig("$savepath/adjustment_by_income.png")
end



function plot_cev_distribution(cev::Array{Float64,5}, dist5::Array{Float64,5}, theta::Float64)
    cev_vec = vec(cev)
    weights = vec(dist5)
    wsum = sum(weights)
    if wsum <= 0
        @warn "Zero weights in plot_cev_distribution"; return
    end
    weights ./= wsum
    histogram(cev_vec, weights=weights, bins=50,
        xlabel="Welfare Gains", ylabel="Density",
        title="Welfare Gains (CEV) — θ = $(round(theta, digits=2))", legend=false)
    savefig("Output/Ergodic/cev_distribution_theta$(round(Int, theta)).png")
end

