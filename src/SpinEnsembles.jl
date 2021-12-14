"""
A julia package that provides necessary functionalities to study the 
quantum decoherence in disordered spin ensembles.
"""
module SpinEnsembles

using LinearAlgebra
using Statistics
import SpecialFunctions: gamma
import ProgressMeter: @showprogress

# datetpye
export SpinEnsemble, SpinCluster

# mathods and functions
export randlocs, randcoefs, 
       fid, averagefid, rabi, betasampling, dephasingtime

# visualization functions
export visualcoupling, visualeffectivebeta, visualensemble, visualfid

# constant
export fid_plot_options

# optional
export dipolarcoef, dipolarcoefs, dipolarlinewidth

include("randloctions.jl")
include("visualization.jl")
include("spindynamics.jl")


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

    function SpinEnsemble(n::Int,dim::Int,z0::AbstractVector{<:Real},
    r::Real, R::Real, shape=:spherical)
        @assert dim in (1,2,3)
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
    SpinEnsemble(ensemble)
    SpinEnsemble(locations)

# Arguments:
- `ensemble`: 
"""
mutable struct SpinCluster
    ensemble::SpinEnsemble
    locations::Matrix{Float64}
    couplings::Vector{Float64}
    linewidth::Float64

    function SpinCluster(ensemble::SpinEnsemble)
        locs=randlocs(ensemble)
        D=dipolarcoefs(locs)
        new(ensemble, locs, D, dipolarlinewidth(D))
    end 

    SpinCluster(locations::Matrix{Float64})=
    (D=dipolarcoefs(locations); new(missing, locations, D, dipolarlinewidth(D)))
end

isdilute(ensemble::SpinEnsemble)=ensemble.rho<1
isdilute(spins::SpinCluster)=spins.ensemble.rho<1

function randlocs(ensemble::SpinEnsemble)
    return randlocs(ensemble.n, ensemble.dim, (ensemble.r,ensemble.R); method=ensemble.shape)
end

randcoefs(ensemble::SpinEnsemble)=dipolarcoefs(randlocs(ensemble.n, ensemble.dim, (ensemble.r,ensemble.R); method=ensemble.shape))


@doc raw"""
    betasampling(D; N)

Give a random sampling of beta, which is a combination of `D_j`,
```math
\beta_p = \sum_j p_j \,D_j,\; p_j=\pm 1
```

# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
- `ensemble`: a `SpinEnsemble`, use this function as the method of the type 
# Options
- `N`: size of Monte-Carlo sampling
"""
function betasampling(D::Vector{<:Real}; N=1::Int)
    n=length(D)
    return [sum(rand([1,-1],n).*D) for i in 1:N]
end

"""
    betasampling(ensemble; N)

Get the Gaussian linewidth of dipolar coupling for the given spin cluster
"""
function betasampling(spins::SpinCluster; N=1::Int)
    return betasampling(spins.couplings, N=N)
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
dipolarlinewidth(spins::SpinCluster)=spins.linewidth # in time computed and stored

@doc raw"""
    dephasingtime(D, n_t=500; scale=1.2)

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
- `n_t`: size of the generated time array 
# Options 
- `scale`: scale factor to extend the T_2 
"""
function dephasingtime(D::Vector{<:Real}, n_t=500::Int; scale=1.2::Real)
    b=dipolarlinewidth(D)
    t=π/b*scale 
    dt=t/n_t
    return 0:dt:t 
end
dephasingtime(spins::SpinCluster, n_t=500::Int; scale=1.2::Real)=dephasingtime(spins.D, n_t; scale=scale)

dephasingtime(D::Vector{<:Real})=π/dipolarlinewidth(D)
dephasingtime(spins::SpinCluster)=π/dipolarlinewidth(spins.D)

"""

"""
function dephasingtime(ensemble::SpinEnsemble)
    if isdilute(ensemble)
        if ensemble.dim==3
            return 6/(ensemble.rho*π^2)*(3*sqrt(3))/8
        elseif ensemble.dim==2
            return 2*(ensemble.rho*π/2*gamma(1/3))^(-3/2)
        elseif ensemble.dim==1
            return -2*(ensemble.rho/sqrt(3)*gamma(-1/3))^(-3)
        else
            throw(DomainError(ensemble.dim,"Invalide ensemble dimension."))
        end
    else
        throw(DomainError("Ensemble is not dilute"))
    end
end

function dephasingtime(ensemble::SpinEnsemble, n_t=500; scale=1.2::Real)
    t=dephasingtime(ensemble)*scale
    dt=t/n_t
    return 0:dt:t 
end

end