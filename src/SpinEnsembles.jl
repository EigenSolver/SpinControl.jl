"""
    SpinEnsembles

A julia package that provides necessary functionalities to study the 
quantum decoherence in disordered spin ensembles.
"""
module SpinEnsembles

using LinearAlgebra
using Statistics
import ProgressMeter: @showprogress

export dipolarcoef, dipolarcoefs, dipolarlinewidth, randomcoefs,
       fid, averagefid, fidsampling, betasampling, decaytime, 
       visualcoupling, visualeffectivebeta, visualensemble, visualfid


include("RandLoctions.jl")
include("Visualization.jl")


"""
    dipolarcoef(r, z0)

Calculate the dipolar interaction strength given by vector `r` and default field at `z0`

# Arguments
- `r::AbstractArray{<:Real}`: point vectors from one loc to another    
- `z0::AbstractArray{<:Real}`: direction of the background magnetic field

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
function dipolarcoef(r::Vector{<:Real},z0::AbstractArray{<:Real})
    # suppose z0 is already normalized
    cosθ=dot(r,z0)/norm(r) #calculate the cos(θ) between the vector and the z axis
    D=0.5*(1-3cosθ^2)/norm(r)^3
    return D
end;


"""
    dipolarcoefs(bath, z0)

Get a list of dipolar coupling strength between the centered spin and bath

# Arguments
- `locs::Matrix{<:Real}`: an array of vector, distance from the central spin to the spins in bath 
- `z0::Vector{<:Real}`: the direction of external field, set to z axis by default
"""
function dipolarcoefs(locs::Matrix{<:Real},z0=[0,0,1.0]::Vector{<:Real})
    normalize!(z0)
    map(x->dipolarcoef(x,z0), eachrow(locs))
end


"""
    randomcoefs(N, dim, a)
    randomcoefs(N, dim, bound, [method])

Randomly distribute a spin bath, generate the dipolar coupling strength between the centered spin and bath, 
totally `N` spins are uniformly distributed in a `dim` dimensional cubic space with lenght `a`

# Arguments
- `N`: number of spin in bath 
- `dim`: dimension
- `a`: scale of ensemble
- `bound::Tuple{Real, Real}`: tuple, indicate bounds of sampling range 
- `method`: constant, `:spherical` for spherical coordinates or `:cubic` for Cartesian coordinates
"""
randomcoefs(N::Int,dim::Int,a=1::Real)=dipolarcoefs(rand_locs(N,dim,a))

function randomcoefs(N::Int, dim::Int, bound::Tuple{Real,Real}; method=:cubic)
    @assert dim>0 && dim<4
    @assert method in (:spherical, :cubic)
    a,b=bound
    @assert a>0 && b>0 

    
    if method==:cubic
        M=randlocscubic(a,b, N=N, dim=dim)
    else
        if dim==3
            M=randlocsspherical(a,b,N=N)
        elseif dim==2
            M=randlocspolar(a,b,N=N)
        else
            M=randlocscubic(a,b,N=N,dim=1)
        end
    end
    dipolarcoefs(M)
end

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
    betasampling(D)

Give a random sampling of beta, which is a combination of `D_j`,
```math
\beta_p = \sum_j p_j \,D_j,\; p_j=\pm 1
```

# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
"""
betasampling(D::Vector{<:Real})=sum(rand([1,-1],length(D)).*D) 

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
    fidsampling(t, D, h; [N])

Get a random sampling of the `f=S_x(t)`, under given transverse magnetic field

```math
f_p(t)=\frac{1}{2}[\cos^2(\omega_p t)+\sin^2(\omega_p t) (n_x^2-n_z^2)]
```

# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 
- `N`: size of Monte-Carlo sampling
"""
function fidsampling(t::AbstractVector{<:Real},D::Vector{<:Real},h::Real;N=1::Int)
    n=length(D)
    f_sum=zeros(length(t)) # sum
    f_var=copy(f_sum) # square sum
    for i in 1:N
        beta_p=sum(rand([1,-1],n).*D) 
        omega_p=sqrt(h^2+beta_p^2)/2
        cos_p=cos.(omega_p*t).^2
        f_p=(cos_p+(cos_p.-1)*(beta_p^2-h^2)/(h^2+beta_p^2))/2
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum).^2/(i*(i-1)) : f_var
    end
    return f_sum/N,f_var/(N-1)
end

function fidsampling(t::Real,D::Vector{<:Real},h::Real;N=1::Int)
    n=length(D)
    f_sum=0.0; f_var=0.0 # sum / square sum
    for i in 1:N
        beta_p=sum(rand([1,-1],n).*D) 
        omega_p=sqrt(h^2+beta_p^2)/2
        cos_p=cos(omega_p*t)^2
        f_p=(cos_p+(cos_p-1)*(beta_p^2-h^2)/(h^2+beta_p^2))/2
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum)^2/(i*(i-1)) : f_var
    end
    return f_sum/N,f_var/(N-1)
end

"""
    decaytime(D, N; len=500, n_sigma=2)

Given a set of coupling strength, determine the average `beta` by sampling.
The average decay time of FID is given by `T=2π/beta`, the `T` is set shorter by 
adding a standard deviation to `beta` given by `n_sigma`
 
# Arguments
- `D::Vector{Real}`: a set of coupling strengths
- `M`: sampling size
- `len`: size of the generated time array 
- `n_sigma`: add standard deviations (σ in Gaussian distribution) to the averaged `beta`
"""
function decaytime(D::Vector{<:Real},M::Int=500;len=500::Int,n_sigma=2::Real)
    sample=abs.([betasampling(D) for i in 1:M])
    omega_m,omega_std=mean(sample),std(sample)
    T=2pi/(omega_m+n_sigma*omega_std)
    return collect(0:T/len:T)
end

@doc raw"""
    averagefid(t, n_ensemble, sampling_D; [options]...)

Calculate the average free induction decay over different ensembles (disorders) 

```math
\bar{f}(t)=\sum_k f(t; \{D_i\}_k)
```

# Arguments
- `t`: the time array for the decay curve.
- `n_ensemble::Integer`: the number of samplings on D set.
- `sampling_D::function`: the function to sample over D
"""
function averagefid(t::AbstractVector{<:Real}, n_ensemble::Int, sampling_D; h=0)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:n_ensemble
        f_d=fid(t, sampling_D())
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    return f_sum/n_ensemble, f_var/(n_ensemble-1)
end

end