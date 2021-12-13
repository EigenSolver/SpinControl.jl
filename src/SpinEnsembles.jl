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
export dipolarcoef, dipolarcoefs, dipolarlinewidth, randlocs, randcoefs,
       fid, averagefid, rabi, betasampling, dephasingtime,
       visualcoupling, visualeffectivebeta, visualensemble, visualfid

# constant
export fid_plot_options


include("randloctions.jl")
include("visualization.jl")
include("dynamics.jl")


"""
    SpinEnsemble(n, dim, z0, r, R)
    SpinEnsemble(n, dim, z0, bound)

# Arguments:
- `n::Int64`: number of spins in the ensemble 
- `dim::Int64`: dimension of the ensemble 
- `z0::Vector{Float64}`: direction of background field 
- `r::Float64`: minimum radius of distribution area 
- `R::Float64`: maximum radius of distribution area
- `rho::Float64`: density of spins
"""
struct SpinEnsemble
    n::Int
    dim::Int
    z0::Vector{Float64}
    
    r::Float64
    R::Float64
    rho::Float64
    shape::Symbol

    # a:: Float64
    # f::Real
    # h::Real
    locs::Matrix{Float64}
    D::Vector{Float64}

    function SpinEnsemble(n::Int,dim::Int,z0::AbstractVector{<:Real},
    r::Real, R::Real, shape=:spherical)
        @assert shape in (:spherical, :cubic)

        if shape==:spherical
            rho=(n/(2*(R-r)),n/(π*(R^2-r^2)),n/(4π/3*(R^3-r^3)))[dim]
        elseif shape==:cubic
            rho=(n/(2*(R-r)),n/(4*(R^2-r^2)),n/(8*(R^3-r^3)))[dim]
        end
        
        return new(n,dim,normalize(z0),Float64(r),Float64(R),rho,shape)
    end

    function SpinEnsemble(n::Int,dim::Int,z0::AbstractVector{<:Real},
        bound::Tuple{Real,Real}, shape=:spherical)
        return SpinEnsemble(n, dim, z0, bound[0], bound[1], shape)
    end
end

"""
    randlocs(n, dim, R)
    randlocs(n, dim, bound; method)

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
function randlocs(n::Int, dim::Int, bound::Tuple{Real,Real}; method=:cubic)
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
    return locs
end

randlocs(n::Int,dim::Int,R=1::Real; method=:cubic)=randlocs(n,dim,(0,R);method=method)

function randlocs(spins::SpinEnsemble)
    return randlocs(spins.n, spins.dim, (spins.r,spins.R); method=spins.shape)
end

function randlocs!(spins::SpinEnsemble)
    spins.locs=randlocs(spins)
end



@doc raw"""
    dipolarcoef(r, z0)

Calculate the dipolar interaction strength given by vector `r` and default field at `z0`

```math
D_{ij}=\frac{1-3\cos^2(\theta_{ij})}{2r_{ij}^3} \gamma^2 \hbar
```

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
    return map(x->dipolarcoef(x,z0), eachrow(locs))
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
    return dipolarcoefs(randlocs(n, dim, bound; method=method))
end

randcoefs(n::Int,dim::Int,R=1::Real; method=:cubic)=randcoefs(n,dim,(0,R);method=method)

randcoefs(spins::SpinEnsemble)=dipolar(randlocs(spins.n, spins.dim, (spins.r,spins.R); method=spins.shape))


function randlocs!(spins::SpinEnsemble)
    spins.locs=randcoefs(spins)
end

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

@doc raw"""
    dephasingtime(D, len=500; scale=1.5)

The FID is Fourier transform of the noise spectrum. 
For a Gaussian noise with linewidth `b`, it's characteristic function is 
```math
f(t)=\int_{\infty}^\infty P(\beta) e^{-i \beta t} d\,\beta=\exp(-b^2 \,t^2/2), \quad P(\beta)=\frac{1}{\sqrt{2\pi b^2}} \exp(-\frac{\beta^2}{2b^2})
```
Thus the decay time is given by 
```math
T_2=\frac{\pi}{b}
```


# Arguments
- `D::Vector{Real}`: a set of coupling strengths
- `len`: size of the generated time array 
# Options 
- `scale`: scale factor to extend the T_2 
"""
function dephasingtime(D::Vector{<:Real}, len=500::Int; scale=1.5::Real)
    b=dipolarlinewidth(D)
    t=π/b*scale 
    dt=t/len
    return 0:dt:t 
end

dephasingtime(D::Vector{<:Real})=π/dipolarlinewidth(D)

dephasingtime(spins::SpinEnsemble)=dephasingtime(spins.D)

end