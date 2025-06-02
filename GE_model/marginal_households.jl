threshold = 1e-3  # Define margin threshold
margin_households = []

for i in 1:sz.nYears
    for j in 1:sz.nFirms
        if abs(alld_adjust[i, j] - alld[i, j]) < threshold && adjustment_indicator[i, j] == 0
            push!(margin_households, (i, j))
        end
    end
end

# Visualization of margin households
scatter(alld[margins[:, 1]], alla[margins[:, 1]], title="Margin Households Pre-Adjustment", label="Margin Households")