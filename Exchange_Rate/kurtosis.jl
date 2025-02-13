function kurtosis(data::Vector{Float64})
    n = length(data)
    mean_data = mean(data)
    std_data = std(data)
    fourth_moment = sum((data .- mean_data).^4) / n
    kurt = fourth_moment / (std_data^4) - 3
    return kurt
end