using LaTeXStrings
using Plots
import LsqFit: curve_fit

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

"""

"""
function visual_FID(t::AbstractArray{<:Real}, FID_curve::AbstractArray{<:Real}; 
    fitting=:false,
    cutoff=1::Int,
    s=2::Real,
    logscale=:false)
    

    f= logscale ? log.(FID_curve) : FID_curve

    scatter(t,f;FID_plot_options..., labels="simulation",
    markershape = :vline,
    markersize = 6,
    markeralpha = 0.9,
    )
    
    if fitting
        f_model(t,p)=-(abs.(t/p[1]).^s).+p[2]
        X=t[cutoff:end]
        Y=log.(FID_curve[cutoff:end])
        p0=[1.0,0.2]
        fit = curve_fit(f_model, X, Y, p0);

        Y_fit=logscale ? f_model(X,fit.param) : exp.(f_model(X,fit.param))

        plot!(X, Y_fit,linestyle=:dash,linewidth=3,label="fitting")
    end 
end