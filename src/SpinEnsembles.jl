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
export fid, averagefid, rabi, averagerabi, betasampling, dephasingtime

# visualization functions
export visualcoupling, visualeffectivebeta, visualensemble, visualfid

# constant
export fid_plot_options

# optional
# export dipolarcoef, dipolarcoefs, dipolarlinewidth

include("randloctions.jl")
include("visualization.jl")
include("spindynamics.jl")
include("properties.jl")


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
    rho::Float64
    dim::Int
    z0::Vector{Float64}
    
    r::Float64
    R::Float64
    shape::Symbol
    T2::Float64

    # a:: Float64
    # f::Real
    # h::Real
    function SpinEnsemble(n::Int,dim::Int,z0::AbstractVector{<:Real},
    r::Real, R::Real, shape=:spherical::Symbol)
        @assert dim in (1,2,3)
        @assert shape in (:spherical, :cubic)
        @assert abs(n-rho*_volume(dim, r, R, shape))<1

        return new(n,dim,normalize(z0),float(r),float(R),rho,shape)
    end

    function SpinEnsemble(n::Int,dim::Int,z0::AbstractVector{<:Real},
    r::Real, R::Real, shape=:spherical::Symbol)
        V=_volume(dim, r, R, shape)
        rho=n/V
        return SpinEnsemble(n,dim,z0,r,R,rho,shape)
    end

    function SpinEnsemble(rho::Float64,dim::Int,z0::AbstractVector{<:Real},
        r::Real, R::Real, shape=:spherical)
        V=_volume(dim, r, R, shape)
        n=floor(Int, rho*V)
        return SpinEnsemble(n,dim,z0,r,R,rho,shape)      
    end

    function SpinEnsemble(n::Int,dim::Int,z0::AbstractVector{<:Real},
        bound::Tuple{Real,Real}, shape=:spherical)
        return SpinEnsemble(n, dim, z0, bound[0], bound[1], shape)
    end
    # function SpinEnsemble(n_or_rho::Union{Int,Float64},dim::Int,z0::AbstractVector{<:Real},
    #     r::Real, R::Real, shape=:spherical::Symbol)
    #         V=_volume(dim, r, R, shape)
    #         if n_or_rho isa Int
    #             n=n_or_rho
    #             rho=n/V
    #         elseif n_or_rho isa Float64
    #             rho=n_or_rho
    #             n=floor(Int, rho*V)
    #         return SpinEnsemble(n,dim,z0,r,R,rho,shape)
    #     end
    
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

randlocs(ensemble::SpinEnsemble)=randlocs(ensemble.n, ensemble.dim, (ensemble.r,ensemble.R); method=ensemble.shape)

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

function averagefid(t::AbstractVector{Float64}, ensemble::SpinEnsemble; M=1000::Int, returnerr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=fid(t, randcoefs(ensemble))
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function averagefid(ensemble::SpinEnsemble;
     M=1000::Int, n_t=200::Int, scale=1.0::Real, returnerr=false)
    T2=dephasingtime(ensemble)*scale
    t=0:T2/n_t:T2 
    return averagefid(t, ensemble; M=M, returnerr=returnerr)
end


function averagerabi(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h::Real;
     M=1000::Int, returnerr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=rabi(t, randcoefs(ensemble), h)
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function averagerabi(ensemble::SpinEnsemble, h::Real; 
    M=1000::Int, n_t=200::Int, scale=1.0::Real, returnerr=false)
    T2=dephasingtime(ensemble)*scale
    t=0:T2/n_t:T2 
    return averagerabi(t, ensemble, h; M=M, returnerr=returnerr)
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

"""
    fid(spins; options)

Calculate the free induction dacay of the given spin cluster
"""
function fid(spins::SpinCluster; n_t=200::Int, scale=1.0::Real)# analytical solution
    T2=dephasingtime(spins)*scale
    t=0:T2/n_t:T2
    return fid(t, spins.couplings)
end

# mote-carlo
function fid(spins::SpinCluster,h::Real; 
    N=100::Int, n_t=200::Int, scale=1.0::Real, geterr=:false)
    T2=dephasingtime(spins)*scale
    t=0:T2/n_t:T2
    fid(t, spins.couplings, h; N=N, geterr=geterr)
end

# analytical solution
function rabi(spins::SpinCluster,h::Real; 
    axis=3::Int, n_t=200::Int, scale=1.0::Real)
    T2=dephasingtime(spins)*scale
    t=0:T2/n_t:T2
    return rabi(t, spins.couplings,h; axis=axis)
end


# mote-carlo
function rabi(spins::SpinCluster, h::Real; 
    N=100::Int, axis=3::Int, n_t=200::Int, scale=1.0::Real, returnerr=false)
    T2=dephasingtime(spins)*scale
    t=0:T2/n_t:T2
    return rabi(t, spins.couplings, h; N=N, axis=axis, returnerr=returnerr)
end

end