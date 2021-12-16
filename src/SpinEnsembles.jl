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
export fid, fid, rabi, rabi, betasampling, coherencetime, randlocs, randcoefs, dipolarlinewidth

include("randloctions.jl")
include("spindynamics.jl")
include("dipolarcouplings.jl")


"""
    SpinEnsemble(n, dim, z0, r, R)
    SpinEnsemble(rho, dim, z0, r, R)

# Arguments:
- `n::Int64`: number of spins in the ensemble 
- `rho::Float64`: density of spins
- `dim::Int64`: dimension of the ensemble 
- `z0::Vector{Float64}`: direction of background field 
- `r::Float64`: minimum radius of distribution area 
- `R::Float64`: maximum radius of distribution area
"""
struct SpinEnsemble
    n::Int
    rho::Float64
    dim::Int
    z0::Vector{Float64}
    r::Float64
    R::Float64
    shape::Symbol

    function SpinEnsemble(n_or_rho::Union{Int,Float64},dim::Int,z0::AbstractVector{<:Real},
        r::Real, R::Real, shape=:spherical::Symbol)
        @assert n_or_rho>0
        @assert dim in (1,2,3)
        @assert shape in (:spherical, :cubic)
        V=_volume(dim, r, R, shape)
        if n_or_rho isa Int
            n=n_or_rho
            rho=n/V
        elseif n_or_rho isa Float64
            rho=n_or_rho
            n=floor(Int, rho*V)
        end
        return new(n,rho,dim,normalize(z0),float(r),float(R),shape)  
    end
end

function _volume(dim::Int, r::Real, R::Real, shape::Symbol)
    if shape==:spherical
        V=((2*(R-r)),(π*(R^2-r^2)),(4π/3*(R^3-r^3)))[dim]
    elseif shape==:cubic
        V=((2*(R-r)),(4*(R^2-r^2)),(8*(R^3-r^3)))[dim]
    end
    return V
end

function volume(ensemble::SpinEnsemble)
    return _volume(ensemble.dim, ensemble.r, ensemble.R, ensemble.shape)
end


# methods & properties
isdilute(ensemble::SpinEnsemble)=ensemble.rho<1
randlocs(ensemble::SpinEnsemble)=randlocs(ensemble.n, ensemble.dim, (ensemble.r,ensemble.R);
 method=ensemble.shape)
randcoefs(ensemble::SpinEnsemble)=dipolarcoefs(randlocs(ensemble))

function coherencetime(ensemble::SpinEnsemble)
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
        throw(DomainError("Ensemble is non-dilute"))
    end
end

mutable struct SpinCluster
    ensemble::SpinEnsemble
    locations::Matrix{Float64}
    couplings::Vector{Float64}
    linewidth::Float64

    function SpinCluster(ensemble::SpinEnsemble)
        locs=randlocs(ensemble)
        D=dipolarcoefs(locs)
        cluster=new(ensemble, locs, D)
        cluster.linewidth=dipolarlinewidth(cluster)
        return cluster
    end 

    function SpinCluster(locations::Matrix{Float64})
        D=dipolarcoefs(locations); 
        cluster=new(missing, locations, D)
        cluster.linewidth=dipolarlinewidth(cluster)
        return cluster
    end
end

isdilute(cluster::SpinCluster)=cluster.ensemble.rho<1

@doc raw"""
    betasampling(cluster; N)

Get the Gaussian linewidth of dipolar coupling for the given spin cluster

    ```math
\beta_p = \sum_j p_j \,D_j,\; p_j=\pm 1
```

# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
- `ensemble`: a `SpinEnsemble`, use this function as the method of the type 
# Options
- `N`: size of Monte-Carlo sampling
"""
function betasampling(cluster::SpinCluster; N=1::Int)
    return [sum(rand([1,-1],n).*cluster.couplings) for i in 1:N]
end



@doc raw"""
    dipolarlinewidth(cluster)

Get the linewidth of D, which follows Gaussian distribution. 

```math
b=\sqrt{\sum_j D_j^2}
```
# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
"""
dipolarlinewidth(cluster::SpinCluster)=sqrt(mapreduce(abs2,+,cluster.couplings)) # in time computed and stored


@doc raw"""
    coherencetime(D, n_t=500; scale=1.0)

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
coherencetime(cluster::SpinCluster)=π/dipolarlinewidth(cluster)

"""
Reroll the points in the cluster
"""
function reroll!(cluster::SpinCluster)
    cluster.locations=randlocs(cluster.ensemble)
    cluster.couplings=dipolarcoefs(cluster.locations)
end

end