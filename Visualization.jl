using LaTeXStrings
using Plots
# using PlotlyJS


"""
Options used to plot a FID line
"""
FID_plot_options=:xformatter=>:scientific,:xlabel=>L"t", :ylabel=>L"$\langle S_x(t) \rangle$",:labels=>:false

"""
Visualize a spin ensemble in 3D dimension
=======
Args:
    spin_locs: (N,3) matrix
Return:
    interactive 3D plot in PlotlyJS
"""
function visual_ensemble(spin_locs::Matrix{Float64})
    plotlyjs()
    scatter(spin_locs[:,1],spin_locs[:,2],spin_locs[:,3],
    xlabel="x",
    ylabel="y",
    zlabel="z",
    markersize=2,
    title="Ensemble Visualization",
    framestyle=:box)
end
