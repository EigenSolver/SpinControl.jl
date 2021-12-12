using LaTeXStrings
using Plots
import LsqFit: curve_fit
import Distributions: Normal, loglikelihood, fit, pdf, params

"""
    plot(...; fid_plot_options...)

Options used to plot a FID line
"""
fid_plot_options=:xformatter=>:scientific,
:xlabel=>L"t", :ylabel=>L"$\langle S_x(t) \rangle$",:labels=>:false

"""
    visualensemble(spin_locs)

Visualize a spin ensemble in 3D dimension, return interactive 3D plot in PlotlyJS

# Arguments
- `spin_locs::Matrix{N,3}`: collections of vectors for `N` spins, in 3 dimension
"""

function visualensemble(spin_locs::Matrix{Float64}; interactive=:true)
    if interactive
        plotlyjs()
    end

    fig=scatter(spin_locs[:,1],spin_locs[:,2],spin_locs[:,3],
    xlabel="x",
    ylabel="y",
    zlabel="z",
    markersize=2,
    title="Ensemble Visualization",
    framestyle=:box)
    
    display(fig);
    gr();
end

"""
    visualfid(t, FID_curve; options...)

Visualize a spin ensemble in 3D dimension, return interactive 3D plot in PlotlyJS

# Arguments
- `t`: time array 
- `FID_curve`: values of FID on the time array `t`

# Options
- `fitting=:false`: do a exponential curve fitting on the FID, set to `:false` by default
- `cutoff=1`: number of points to drop at the beginning of fitting 
- `s=2`: exponential power of the fitting model 
- `logscale=:false`: display the plot in logscale, set to `false` by default 
"""
function visualfid(t::AbstractArray{<:Real}, FID_curve::AbstractArray{<:Real}; 
    fitting=:false,
    cutoff=1::Int,
    s=2::Real,
    logscale=:false)
    

    f= logscale ? log.(FID_curve) : FID_curve

    fig = scatter(t,f;fid_plot_options..., labels="simulation",
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



"""
    visualcoupling(sample; options...)

Visualize an array of coupling strengths using histogram

# Arguments
- `sample`: an array of coupling strengths

# Options
- `bin_set`: bins in the histogram
- `logscale`: whether to display the histogram in logscale
"""
function visualcoupling(sample::AbstractArray{Float64}; 
    bin_set=:none, logscale=:false)
    
    if bin_set==:none
        bin_set=range(-20, 20, length = 60)
    end 

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


"""
    visualeffectivebeta(sample; options...)

Visualize an array of effective magnetic field 

# Arguments
- `sample`: an array of beta

# Options
- `bin_set`: bins in the histogram
- `use_abs`: whether to use the absolute value of `sample` in the histogram
- `fitting`: fit the distribution with Gaussian distribution
"""
function visualeffectivebeta(sample::AbstractArray{Float64};
    bin_set=-50:2:50, use_abs=:false, fitting=:true)
    
    if use_abs 
        sample=abs.(sample)
    end

    fig=histogram(sample, bins = bin_set, xlabel=L"\beta_p", ylabel=L"\Delta P", norm=false,label="Beta Sample")
    println("Distribution:")
    println("max: ",maximum(sample)," min: ",minimum(sample))
    println("avg: ", mean(sample)," std: ",std(sample))
    
    if fitting
        normal_est=fit(Normal{Float64},sample)
        mu,sigma=params(normal_est)
        # bins=LinRange(mu-4sigma,mu+4sigma,bin_num*2+1)
        plot!(bins,map(x->pdf(normal_est,x),bin_set),linecolor=:red,linestyle=:dash,label="Normal Distribution")
        println("Estimation:")
        println("mu: ", mu, " sigma: ", sigma)
        println("likelihood: ", ℯ^loglikelihood(normal_est,sample))
    end

    display(fig)
end;