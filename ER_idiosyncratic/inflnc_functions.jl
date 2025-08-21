# Mean influence function
function mean_inflnc(z);
    meaninflnc = z .- mean(z);
    return meaninflnc;
end;

# Variance influence function
function var_inflnc(z)
    vz = var(z)
    varinflnc = (z .- mean(z)).^2.0 .- vz;
    return varinflnc;
end

# Covariance influence function
function cov_inflnc(xy)
    covxy =  mean(x.*y) - mean(x)*mean(y)
    covinflnc = (x .- mean(x)).*(y .- mean(y)) .- covxy;
    return covinflnc;
end

# OLS influence function 
function ols_inflnc(y,x);
    n = size(y,1);
    bhat = do_ols(y,x); 
    uhat = y - x*bhat; 
    olsinflnc = (inv((x'*x)./n) * ((x.*uhat)'))';
    return olsinflnc;
end; 

#OLS function
function do_ols(y,x)
    bhat = inv(x'*x)*x'*y;
    return bhat;
end

#Cluster function
function cluster(f,g)
    #Inputs:    f is a matrix with dimension of # of obs by # of influence functions
    #           g is the clustering variable (state, year, company, etc.)
    gi = unique(g)
    k = size(f,2);
    vmx = zeros(k,k);

    for ig in eachindex(gi); #for all of the unique identifiers 
        phii = f[g.==gi[ig],:];
        phii = sum(phii,dims=1); #println(phii'*phii)
        vmx = vmx + phii'*phii;
    end

    vmx = vmx./(size(f,1)^2.0);
    return vmx;

end; 

function doublecluster(f,g1,g2)
    covmx1 = cluster(f,g1)
    covmx2 = cluster(f,g2)
    covmx3 = f'*f./(size(f,1)^2.0); 
    covmx = covmx1 .+ covmx2 .- covmx3;
    return covmx;
end;

function ratio_inflnc(x::Vector{Float64}, y::Vector{Float64})
    mx = mean(x)
    my = mean(y)
    infl = (x .- mx)./my .- (mx./(my^2)).*(y .- my)
    return infl
end

function ratio_inflnc(x::Vector{T}) where T<:Real
    R = mean(x .> 0)
    return (x .> 0) .- R
end