"""
    SpinEnsembles

A julia package that provides necessary functionalities to study the 
quantum decoherence in disordered spin ensembles.
"""
module SpinEnsembles

using LinearAlgebra
using Statistics
import ProgressMeter: @showprogress

# functions
export dipolarcoef, dipolarcoefs, dipolarlinewidth, randcoefs,
       fid, averagefid, rabi, betasampling, decaytime,
       visualcoupling, visualeffectivebeta, visualensemble, visualfid

# constant
export fid_plot_options


include("RandLoctions.jl")
include("Visualization.jl")


"""
    SpinEnsemble
 
"""
struct SpinEnsemble
    n::Int  
    dim::Int 
    h::Real
    z0::AbstractVector{<:Real}
    
    r::Real
    R::Real
    a:: Real
    rho::Real
    # f::Real

    locs::AbstractMatrix{<:Real}
    D::AbstractVector{<:Real}

    # islatticed::Bool
    # isdilute::Bool
end

"""
    dipolarcoef(r, z0)

Calculate the dipolar interaction strength given by vector `r` and default field at `z0`

# Arguments
- `r`: point vectors from one loc to another    
- `z0`: direction of the background magnetic field

# Examples
```jldoctest
julia> v0 = [2, 2, 1]; z0=[0, 0, 1];
julia> dipolarcoef(v0,z0)
0.01234567901234568
julia> v0 = [1, 0, 1]; z0=[0, 0, 1];
julia> dipolarcoef(v0,z0)
-0.08838834764831834
```    
"""
function dipolarcoef(r::AbstractVector{<:Real},z0::AbstractVector{<:Real})
    # suppose z0 is already normalized
    cosθ=dot(r,z0)/norm(r) #calculate the cos(θ) between the vector and the z axis
    D=0.5*(1-3cosθ^2)/norm(r)^3
    return D
end


"""
    dipolarcoefs(locs, z0)

Get a list of dipolar coupling strength between the centered spin and bath

# Arguments
- `locs`: an array of vector, distance from the central spin to the spins in bath 
- `z0`: the direction of external field, set to z axis by default
"""
function dipolarcoefs(locs::Matrix{<:Real},z0=[0,0,1.0]::AbstractVector{<:Real})
    normalize!(z0)
    map(x->dipolarcoef(x,z0), eachrow(locs))
end


"""
    randcoefs(n, dim, R)
    randcoefs(n, dim, bound; method)

Randomly distribute a spin bath, generate the dipolar coupling strength between the centered spin and bath, 
totally `n` spins are uniformly distributed in a `dim` dimensional cubic space with lenght `a`

# Arguments
- `n`: number of spins in ensemble
- `dim`: dimension
- `R`: scale of ensemble
- `bound::Tuple{Real, Real}`: tuple, indicate bounds of sampling range 

# Options
- `method`: constant, `:spherical` for spherical coordinates or `:cubic` for Cartesian coordinates
"""
function randcoefs(n::Int, dim::Int, bound::Tuple{Real,Real}; method=:cubic)
    @assert dim>0 && dim<4
    @assert method in (:spherical, :cubic)
    a,b=bound
    @assert a>=0 && b>=0
    a,b= a<b ? (a,b) : (b,a)

    if method==:cubic
        locs=randcartesianlocs(n,a,b, dim=dim)
    else
        if dim==3
            locs=randsphericallocs(n,a,b)
        elseif dim==2
            locs=randpolarlocs(n,a,b)
        else
            locs=randcartesianlocs(n,a,b, dim=1)
        end
    end
    dipolarcoefs(locs)
end

randcoefs(n::Int,dim::Int,R=1::Real; method=:cubic)=randcoefs(n,dim,(0,R);method=method)

@doc raw"""
    betasampling(D)

Give a random sampling of beta, which is a combination of `D_j`,
```math
\beta_p = \sum_j p_j \,D_j,\; p_j=\pm 1
```

# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
# Options
- `N`: size of Monte-Carlo sampling
"""
function betasampling(D::Vector{<:Real}; N=1::Int)
    n=length(D)
    return [sum(rand([1,-1],n).*D) for i in 1:N]
end

@doc raw"""
    dipolarlinewidth(D)

Get the linewidth of D, which follows Gaussian distribution. 

```math
b=\sqrt{\sum_j D_j^2}
```
# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
"""
dipolarlinewidth(D::Vector{<:Real})=sqrt(mapreduce(abs2,+,D))

@doc raw"""
    fid(t, D)

Given a set of dipolar coupling constant, 
calculate et the free induction decay (FID) value at given time `t`, or the FID decay curve for a time array `t`

```math
f(t)=\frac{1}{2}\prod_j \cos(D_jt)
```
# Arguments
- `t`: time, float or array of float
- `D`: coupling strength of dipolar interactions, a array of floats
"""
fid(t::Real,D::Vector{<:Real})=mapreduce(cos,*,D*t)/2;
fid(t::AbstractVector{<:Real},D::Vector{<:Real})=map(x->fid(x,D),t);

@doc raw"""
    fid(t, D, h; N)

Simulate FID signal with given transverse magnetic field `h`, using Monte-Carlo sampling

```math
f_p(t)=\frac{1}{2}[\cos^2(\omega_p t)+\sin^2(\omega_p t) (n_x^2-n_z^2)]
```

```math
\bar{f}(t)=\frac{1}{N}\sum_{p=1}^N f_p(t)
```
# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `N`: number of Monte-Carlo sampling
- `geterr::Bool`: wetehr to return the error of the monte-sampling 
"""
function fid(t::AbstractVector{<:Real},D::Vector{<:Real},h::Real; N=100::Int, geterr=:false)
    n=length(D)
    f_sum=zeros(length(t)) # sum
    f_var=copy(f_sum) # square sum
    for i in 1:N
        β_p=sum(rand([1,-1],n).*D) 
        ω_p=sqrt(h^2+β_p^2)/2
        cos2_p=cos.(ω_p*t).^2
        f_p=(cos2_p+(cos2_p.-1)*(β_p^2-h^2)/(h^2+β_p^2))/2
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum).^2/(i*(i-1)) : f_var
    end
    if geterr
        return (f_sum/N,f_var/(N-1))
    else
        return f_sum/N
    end
end


@doc raw"""
    averagefid(t, M, sampling_D; [options]...)

Calculate the average free induction decay over different ensembles (disorders) 

```math
\bar{f}(t)=\sum_k f(t; \{D_i\}_k)
```

# Arguments
- `t`: the time array for the decay curve.
- `M`: number of ensembles  
- `sampling_D`: the function to sample over D
"""
function averagefid(t::AbstractVector{<:Real}, M::Int, sampling_D::Function)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=fid(t, sampling_D())
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    return f_sum/M, f_var/(M-1)
end

"""
    decaytime(D, N; len=500, n_sigma=2)

Given a set of coupling strength, determine the average `beta` by sampling.
The average decay time of FID is given by `T=2π/beta`, the `T` is set shorter by 
adding a standard deviation to `beta` given by `n_sigma`

# Arguments
- `D::Vector{Real}`: a set of coupling strengths
- `N`: sampling size
- `len`: size of the generated time array 
- `n_sigma`: add standard deviations (σ in Gaussian distribution) to the averaged `beta`
"""
function decaytime(D::Vector{<:Real},N::Int=500; len=500::Int, n_sigma=2::Real)
    sample=abs.(betasampling(D;N=N))
    omega_m,omega_std=mean(sample),std(sample)
    T=2pi/(omega_m+n_sigma*omega_std)
    return collect(0:T/len:T)
end

@doc raw"""
    rabi(t, D, h; N, options...)

Get a random sampling of, under given transverse magnetic field

```math
\begin{aligned}
z_p(t)&=+\frac{1}{2} \left[n_z^2+n_x^2 \cos (\Omega_p t)\right]\\
y_p(t)&=-\frac{1}{2} n_x \sin (\Omega_p t)\\
x_p(t)&=+\frac{1}{2} \left[n_x n_z-n_x n_z \cos (\Omega_p t)\right]
\end{aligned}
```

```math
G(t)=\frac{1}{N}\sum_{p=1}^N g_p(t),\; g=x,y,z
```

# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `N`: size of Monte-Carlo sampling, default at 100
- `axis::Int`: , 1,2,3, representing x,y,z axis, set to 3 by default 
- `returnerr::Bool`: wether to return the variance in sampling
"""
function rabi(t::AbstractVector{<:Real}, D::Vector{<:Real}, h::Real; N=100::Int, axis=3::Int, returnerr=:false)
    @assert axis in (1,2,3)
    n=length(D)
    f_sum=zeros(length(t)) # sum
    f_var=copy(f_sum) # square sum
    f_sampling=(_rabix,_rabiy,_rabiz)[axis]

    for i in 1:N
        β_p=sum(rand([1,-1],n).*D) 
        f_p=f_sampling(t,β_p,h)
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return (f_sum/N,f_var/(N-1))
    else
        return f_sum/N
    end
end

function _rabiz(t::AbstractVector{<:Real},β::Real,h::Real)
    ω=sqrt(h^2+β^2)/2
    cos2=cos.(ω*t).^2
    return (cos2-(cos2.-1)*(β^2-h^2)/(h^2+β^2))/2
end

function _rabiy(t::AbstractVector{<:Real},β::Real,h::Real)
    Ω=sqrt(h^2+β^2)
    return -sin.(Ω*t)/2*h/Ω
end

function _rabix(t::AbstractVector{<:Real},β::Real,h::Real)
    ω=sqrt(h^2+β^2)/2
    sin2=sin.(ω*t).^2
    return sin2*(β*h)/(h^2+β^2)
end


@doc raw"""
    rabi(t, h, b; [axis])

Get analytical solution of Rabi oscillation under Gaussian noise with linewidth `b`

```math
M_z(t)=A \cos(h \,t+\varphi), \quad M_y(t)=-A \sin(h \,t+\varphi)
```
```math
A=\frac{1}{2}\left(1+b^4t^2/h^2\right)^{1/4},\quad \varphi=\frac{1}{2}\arctan(b^2 t/h)
```

# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `axis::Int`: , 1,2,3, representing x,y,z axis, set to 3 by default 
""" 
function rabi(t::AbstractArray, h::Real, b::Real; axis=3)
    @assert axis in (1,2,3)
    A=(1.0.+b^4*t.^2/h^2).^(-1/4)/2
    φ=atan.(b^2*t/h)/2
    if axis==3
        return A.*cos.(h*t+φ)
    elseif axis==2
        return -A.*sin.(h*t+φ)
    end
end

end