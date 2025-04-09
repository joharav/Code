using Plots
gr() # switch to the GR backend
default(fontfamily = "Computer Modern")  # Looks like LaTeX

# Define function to calculate average moment across variations for a given parameter
function avg_moments_for_param(allmoms, param_index, nmom)
    averages = [mean(allmoms[:, param_index, imom]) for imom in 1:nmom]
    return averages
end

# Assuming you have 4 parameters and 4 moments as per your description
nparam = sz.nop
nmom = sz.nmom

# Initialize matrix to hold average moment values for the heatmap
heatmap_matrix = zeros(nmom, nparam)

# Compile data for heatmap
for iparam in 1:nparam
    heatmap_matrix[:, iparam] .= avg_moments_for_param(allmoms, iparam, nmom)
end

# Now we can create the heatmap
heatmap_plot = heatmap(
    pname, # parameter names
    momname, # moments names
    heatmap_matrix,
    cmap = :blues, # Specify the blue color gradient
    aspect_ratio=:equal,
    xlabel="Parameters",
    ylabel="Moments",
    title="Heatmap of Parameters vs. Moments"
)

# Display the plot
display(heatmap_plot)

# Save the heatmap to a file if needed
savefig(heatmap_plot, "Output/Comparative/Parameters_vs_Moments_Heatmap.png")