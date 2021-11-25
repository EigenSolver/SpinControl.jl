using LaTeXStrings
using Plots
import LsqFit: curve_fit

# plotlyjs()

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

    fig = scatter(t,f;FID_plot_options..., labels="simulation",
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
        println("T_2 = ",fit.param[1])
        println("Residual: ", mean(fit.resid))

        plot!(X, Y_fit,linestyle=:dash,linewidth=3,label=L"fitting $e^{-(t/T_2)^{d/3}}$")
    end 
    display(fig)
end


function visual_coupling(sample::AbstractArray{Float64}; 
    bound=(-20,20), bin_num=60::Int, logscale=:false, fitting=:false)
    
    bin_set=range(bound[1], bound[2], length = bin_num)
    if logscale
        fig=histogram(sample, bins = bin_set, xlabel=L"D_j", ylabel=L"\Delta P", 
        yaxis =:log10,
        norm=true, labels=:false)
    else
        fig=histogram(sample, bins = bin_set, xlabel=L"D_j", ylabel=L"\Delta P", 
        norm=true, labels=:false)
    end
    println("max: ",maximum(sample)," min: ",minimum(sample))
    println("avg: ", mean(sample)," std: ",std(sample))
    display(fig)
end;


function visual_effective_beta(sample::AbstractArray{Float64}; bin_num=60::Int, use_abs=:false, fitting=:true)
    if use_abs 
        sample=abs.(sample)
    end

    fig=histogram(sample, bins = bin_num, xlabel=L"\beta_p", ylabel=L"\Delta P", norm=true,label="Beta Sample")
    println("max: ",maximum(sample)," min: ",minimum(sample))
    println("avg: ", mean(sample),"std: ",std(sample))
    
    if fitting
        normal_est=fit(Normal{Float64},sample)
        mu,sigma=params(normal_est)
        bins=LinRange(mu-4sigma,mu+4sigma,bin_num*2+1)
        plot!(bins,pdf(normal_est,bins),linecolor=:red,linestyle=:dash,label="Normal Distribution")
        println("likelihood: ", ℯ^loglikelihood(normal_est,sample))
    end

    display(fig)
end;