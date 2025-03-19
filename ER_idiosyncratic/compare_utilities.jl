# Function to read utility values from a file into an array
function read_utility_from_file(filename::String, dims)
    util_array = Array{Float64}(undef, dims...)
    i = 1
    open(filename, "r") do io
        for line in eachline(io)
            values = parse.(Float64, split(line))
            for v in values
                util_array[i] = v
                i += 1
            end
        end
    end
    return util_array
end

# Example function comparison between two utility arrays
function compare_utilities(util_noadj::Array{Float64}, util_adj::Array{Float64})
    total_cases = length(util_noadj)
    greater_count = 0

    for idx in eachindex(util_noadj)
        if util_noadj[idx] > util_adj[idx]
            greater_count += 1
        end
    end

    share = greater_count / total_cases * 100
    println("Share of cases where utility_noadjust > utility_adjust: ", share, "%")
end

# Reading utility arrays from files
dims = (sz.ne, sz.na, sz.nd, sz.npa, sz.npd)  # Adjust to appropriate dimensions

util_noadjust = read_utility_from_file("Output/Policy/U_NoAdjust.txt", dims)
util_adjust = read_utility_from_file("Output/Policy/U_Adjust.txt", dims)

# Comparing the utilities
compare_utilities(util_noadjust, util_adjust)