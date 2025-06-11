function simmodel_mit(answ::NamedTuple, e_shock::Float64, T::Int)
    v, pol, grids = answ.v, answ.pol, answ.g

    # Steady state values
    e_ss = 1.0
    z_ss = 1.0
    a0 = 10.0    # or mean from ergodic dist
    d0 = 20.0    # same

    e_path = vcat([e_shock], ones(T - 1))   # Shock only at t=0
    z_path = ones(T)                        # Assume constant productivity

    # Initialize containers
    a = zeros(T)
    d = zeros(T)
    c = zeros(T)
    vval = zeros(T)

    # Initial values
    a[1] = a0
    d[1] = d0
    vval[1] = interpol(z_ss, e_path[1], a0, d0, grids, v)
    c[1] = interpol(z_ss, e_path[1], a0, d0, grids, pol.c)

    for t in 2:T
        a[t] = interpol(z_path[t-1], e_path[t-1], a[t-1], d[t-1], grids, pol.a)
        d[t] = interpol(z_path[t-1], e_path[t-1], a[t-1], d[t-1], grids, pol.d)
        c[t] = interpol(z_path[t-1], e_path[t-1], a[t-1], d[t-1], grids, pol.c)
        vval[t] = interpol(z_path[t-1], e_path[t-1], a[t-1], d[t-1], grids, v)
    end

    return (a=a, d=d, c=c, v=vval, e=e_path, z=z_path)
end


function plot_mit(simul::NamedTuple)
    T = length(simul.c)
    quarters = 0:T-1
    plot(quarters, 100 .* (simul.c ./ simul.c[end] .- 1), label="Consumption", xlabel="Quarter", ylabel="% Deviation from SS")
    plot!(quarters, 100 .* (simul.a ./ simul.a[end] .- 1), label="Assets")
    plot!(quarters, 100 .* (simul.d ./ simul.d[end] .- 1), label="Durables")
end

function compute_elasticity(path_shock::Vector{Float64}, path_noshock::Vector{Float64}, eps::Float64)
    elasticity = (log.(path_shock) .- log.(path_noshock)) ./ eps
    return elasticity
end
