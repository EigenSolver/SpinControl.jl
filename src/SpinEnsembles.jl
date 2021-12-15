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
export fid, averagefid, rabi, averagerabi, betasampling, dephasingtime, randlocs, randcoefs

# visualization functions
export visualcoupling, visualeffectivebeta, visualensemble, visualfid

# constant
export fid_plot_options

# optional
export dipolarcoef, dipolarcoefs, dipolarlinewidth

include("randloctions.jl")
include("visualization.jl")
include("spindynamics.jl")
include("properties.jl")


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
        # @assert abs(n-rho*_volume(dim, r, R, shape))<1

        return new(n,rho,dim,normalize(z0),float(r),float(R),shape)  
    end
    
end

function _volume(dim::Int, r::Real, R::Real, shape::Symbol)
    if shape==:spherical
        V=((π*(R^2-r^2)),(4π/3*(R^3-r^3)))[dim]
    elseif shape==:cubic
        V=((2*(R-r)),(4*(R^2-r^2)),(8*(R^3-r^3)))[dim]
    end
    return V
end

function volume(ensemble::SpinEnsemble)
    return _volume(ensemble.dim, ensemble.r, ensemble.R, ensemble.shape)
end



isdilute(ensemble::SpinEnsemble)=ensemble.rho<1

randlocs(e::SpinEnsemble)=randlocs(e.n, e.dim, (e.r,e.R); method=e.shape)

randcoefs(ensemble::SpinEnsemble)=dipolarcoefs(randlocs(ensemble))


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

function dephasingtime(ensemble::SpinEnsemble, n_t::Int; scale=1.0::Real)
    T2=dephasingtime(ensemble)*scale
    return 0:T2/n_t:T2 
end

function averagefid(t::AbstractVector{Float64}, ensemble::SpinEnsemble; M=1000::Int, geterr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=fid(t, randcoefs(ensemble))
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if geterr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function averagefid(ensemble::SpinEnsemble;
    M=1000::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=dephasingtime(ensemble, n_t; scale=scale)
    return averagefid(t, ensemble; M=M, geterr=geterr)
end

function averagerabi(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h::Real;
    M=1000::Int, geterr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=rabi(t, randcoefs(ensemble), h)
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if geterr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function averagerabi(ensemble::SpinEnsemble, h::Real; 
    M=1000::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=dephasingtime(ensemble, n_t; scale=scale)
    return averagerabi(t, ensemble, h; M=M, geterr=geterr)
end

mutable struct SpinCluster
    n::Int
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

isdilute(spins::SpinCluster)=spins.ensemble.rho<1

"""
    betasampling(ensemble; N)

Get the Gaussian linewidth of dipolar coupling for the given spin cluster
"""
function betasampling(spins::SpinCluster; N=1::Int)
    return [sum(rand([1,-1],n).*spins.couplings) for i in 1:N]
end

dipolarlinewidth(spins::SpinCluster)=spins.linewidth # in time computed and stored

dephasingtime(spins::SpinCluster)=π/dipolarlinewidth(spins.D)

function dephasingtime(spins::SpinCluster, n_t::Int; scale=1.0::Real)
    T2=dephasingtime(spins)*scale
    return 0:T2/n_t:T2 
end


"""
    fid(spins; options)

Calculate the free induction dacay of the given spin cluster
"""
function fid(t::AbstractVector{Float64}, spins::SpinCluster)# analytical solution
    D=spins.couplings
    return [mapreduce(cos,*,D*τ)/2 for τ in t]
end

function fid(spins::SpinCluster; n_t=200::Int, scale=1.0::Real)# analytical solution
    t=dephasingtime(spins,n_t; scale=scale)
    return fid(t,spins)
end

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
function fid(t::AbstractVector{Float64}, spins::SpinCluster,h::Real; 
    N=100::Int, geterr=false)
    D=spins.couplings
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

function fid(spins::SpinCluster,h::Real; 
    N=100::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=dephasingtime(spins,n_t; scale=scale)
    fid(t, spins, h; N=N, geterr=geterr)
end


# analytical solution (approximate at h>>b)
function rabi(t::AbstractVector{Float64}, spins::SpinCluster, h::Real; 
    axis=3::Int)
    @assert axis in (1,2,3)

    b=spins.linewidth
    A=(1.0.+b^4*t.^2/h^2).^(-1/4)/2
    φ=atan.(b^2*t/h)/2
    if axis==3
        return A.*cos.(h*t+φ)
    elseif axis==2
        return -A.*sin.(h*t+φ)
    end
end


# mote-carlo
function rabi(spins::SpinCluster, h::Real; 
    N=100::Int, axis=3::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=dephasingtime(spins,n_t; scale=scale)
    return rabi(t, spins, h; N=N, axis=axis, geterr=geterr)
end

end