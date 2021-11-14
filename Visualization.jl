using LaTeXStrings
using Plots


"""
Options used to plot a FID line
"""
FID_plot_options=:xformatter=>:scientific,:xlabel=>L"t", :ylabel=>L"$\langle S_x(t) \rangle$",:labels=>:false

"""

"""
function visual_ensemble(spin_locs::Matrix{Float64})
    plotlyjs()
    scatter(spin_locs[:,1],spin_locs[:,2],spin_locs[:,3],
    xlabel="x")
end
