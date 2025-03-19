function tauchen(mew::Float64, sigma::Float64, rho::Float64, znum::Int, q::Float64)

    # PACKAGES
    # Distributions

    # INPUT
    # mew the intercept of the AR1 process
    # rho is the AR1 coefficient
    # sigma is the standard deviation of the residual of the AR1 process
    # q is the number of standard deviations from the mean. It really is a float
    # znum is the length of the discretized output vector
    
    # OUTPUT
    # z is the discretized output vector
    # trans is the transition matrix 

    zstar=mew/(1.0-rho); #expected value of z
    sigmaz=sigma/sqrt(1.0-rho^2); #stddev of z

    z = zstar .+ collect(range(-q*sigmaz,stop=q*sigmaz,length=znum));

    trans = zeros(znum,znum)
    w = (z[2] - z[1]);  #Note that all the points are equidistant by construction.
    Threads.@threads for iz = 1:znum
        Threads.@threads for izz = 1:znum
            binhi = (z[izz] - rho * z[iz] + w / 2.0) / sigmaz
            binlo = (z[izz] - rho * z[iz] - w / 2.0) / sigmaz
            if izz == 1
                trans[iz, izz] = cdf(Normal(),binhi)
            elseif izz == znum
                trans[iz, izz] = 1.0 - cdf(Normal(),binlo)
            else
                trans[iz, izz] = cdf(Normal(),binhi) - cdf(Normal(),binlo)
            end
        end
    end

    return z::Vector{Float64}, trans::Matrix{Float64};
end
