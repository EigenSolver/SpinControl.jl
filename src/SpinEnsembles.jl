module SpinEnsembles

using LinearAlgebra
using Statistics
import ProgressMeter: @showprogress

export dipolar_coef, bath_dipolar_coefs, rand_bath_dipolar_coefs,
       f_sampling, beta_sampling, t_adaptive, 
       ensemble_FID, ensemble_average_FID,
       transverse_threshold, 
       visual_coupling, visual_effective_beta, visual_ensemble, visual_FID


include("RandLoctions.jl")
include("Visualization.jl")


"""
    bath_vectors(vec0, vec_bath)

Given the locations of central spin and its bath spin, return the vector set from the central spin to bath.

# Arguments
- `loc0::Vector{<:Real}`, location of the central spin
- `loc_bath::Matrix{<:Real}`: matrix, location set of the spin bath

# Examples
```jldoctest
julia> v0 = [1, 1, 1]; bath= [1.5 1.5 1.0; 2.1 2.2 2.5];
julia> bath_vectors(v0, bath)
2-element Vector{Vector{Float64}}:
 [0.5, 0.5, 0.0]
 [1.1, 1.4, 1.5]
```    
"""
bath_vectors(loc0::Vector{<:Real},loc_bath::Matrix{<:Real})=map(x->x-loc0,eachrow(loc_bath))

"""
    dipolar_coef(r, z0)

Calculate the dipolar interaction strength given by vector `r` and default field at `z0`

# Arguments
- `r::AbstractArray{<:Real}`: point vectors from one loc to another    
- `z0::AbstractArray{<:Real}`: direction of the background magnetic field

# Examples
```jldoctest
julia> v0 = [2, 2, 1]; z0=[0, 0, 1];
julia> dipolar_coef(v0,z0)
0.01234567901234568
julia> v0 = [1, 0, 1]; z0=[0, 0, 1];
julia> dipolar_coef(v0,z0)
-0.08838834764831834
```    
"""
function dipolar_coef(r::AbstractArray{<:Real},z0::AbstractArray{<:Real})
    # suppose z0 is already normalized
    cosθ=dot(r,z0)/norm(r) #calculate the cos(θ) between the vector and the z axis
    D=0.5*(1-3cosθ^2)/norm(r)^3
    return D
end;


"""
    bath_dipolar_coefs(bath, z0)

Get a list of dipolar coupling strength between the centered spin and bath

# Arguments
- `vec_bath::Matrix{<:Real}`: an array of vector, distance from the central spin to the spins in bath 
- `z0::Vector{<:Real}`: the direction of external field, set to z axis by default
"""
function bath_dipolar_coefs(vec_bath::Matrix{<:Real},z0=[0,0,1.0]::Vector{<:Real})
    normalize!(z0)
    map(x->dipolar_coef(x,z0), eachrow(vec_bath))
end


"""
    rand_bath_dipolar_coefs(N, dim, a)
    rand_bath_dipolar_coefs(N, dim, bound, [method])

Randomly distribute a spin bath, generate the dipolar coupling strength between the centered spin and bath, 
totally `N` spins are uniformly distributed in a `dim` dimensional cubic space with lenght `a`

# Arguments
- `N`: number of spin in bath 
- `dim`: dimension
- `a`: scale of ensemble
- `bound::Tuple{Real, Real}`: tuple, indicate bounds of sampling range 
- `method`: constant, `:spherical` for spherical coordinates or `:cubic` for Cartesian coordinates
"""
rand_bath_dipolar_coefs(N::Int,dim::Int,a=1::Real)=bath_dipolar_coefs(rand_locs(N,dim,a))

function rand_bath_dipolar_coefs(N::Int, dim::Int, bound::Tuple{Real,Real}; method=:cubic)
    @assert dim>0 && dim<4
    @assert method in (:spherical, :cubic)
    a,b=bound
    @assert a>0 && b>0 

    
    if method==:cubic
        M=rand_locs_cubic(a,b, N=N, dim=dim)
    else
        if dim==3
            M=rand_locs_spherical(a,b,N=N)
        elseif dim==2
            M=rand_locs_polar(a,b,N=N)
        else
            M=rand_locs_cubic(a,b,N=N,dim=1)
        end
    end
    bath_dipolar_coefs(M)
end

"""
# Arguments
    t: time, float or array of float
    D: coupling strength of dipolar interactions, usually a array of floats
Returns:
    ensemble free induction decay curve at given time t
"""
ensemble_FID(t::Real,D::Vector{<:Real})=mapreduce(cos,*,D*t)/2;
ensemble_FID(t::AbstractVector{<:Real},D::Vector{<:Real})=map(x->ensemble_FID(x,D),t);


"""
Find the amplitude of tranverse magnetic field being stronger than given threshold of the coupling strength
 
# Arguments
    f: density of the spins
    a: 1d scale of the ensemble
    p: probability threshold
"""
function transverse_threshold(p::Real, f::Real, d::Int, a::Real)
    N=floor(Int,a^d*f)
    D=rand_bath_dipolar_coefs(N,d,a)
    quantile(D, p)
end;


"""
Give a random sampling of beta, which is a combination of D
 
# Arguments
    D_set: a set of the coupling strengths
"""
beta_sampling(D_set::Vector{<:Real})=sum(rand([1,-1],length(D_set)).*D_set) 

"""
Get a random sampling of the f=Sx(t), under given transverse magnetic field
 
# Arguments
    t: discrete array marking time 
    D_set: a set of the coupling strengths
    h: strength of transverse field 
    N: size of Monte-Carlo sampling
"""

function f_sampling(t::AbstractVector{<:Real},D_set::Vector{<:Real},h::Real;N=1::Int)
    n=length(D_set)
    f_sum=zeros(length(t)) # sum
    f_var=copy(f_sum) # square sum
    for i in 1:N
        beta_p=sum(rand([1,-1],n).*D_set) 
        omega_p=sqrt(h^2+beta_p^2)/2
        cos_p=cos.(omega_p*t).^2
        f_p=(cos_p+(cos_p.-1)*(beta_p^2-h^2)/(h^2+beta_p^2))/2
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum).^2/(i*(i-1)) : f_var
    end
    return f_sum/N,f_var/(N-1)
end

function f_sampling(t::Real,D_set::Vector{<:Real},h::Real;N=1::Int)
    n=length(D_set)
    f_sum=0.0; f_var=0.0 # sum / square sum
    for i in 1:N
        beta_p=sum(rand([1,-1],n).*D_set) 
        omega_p=sqrt(h^2+beta_p^2)/2
        cos_p=cos(omega_p*t)^2
        f_p=(cos_p+(cos_p-1)*(beta_p^2-h^2)/(h^2+beta_p^2))/2
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum)^2/(i*(i-1)) : f_var
    end
    return f_sum/N,f_var/(N-1)
end

"""
Given a set of coupling strength, determine the maximum time scale and minimum time scale required for the problem
 
# Arguments
    D: a set of coupling strengths
    M: sampling size
    len: size of the generated time array 
Return:
    an array of time points
"""
function t_adaptive(D::Vector{<:Real},M::Int=500;len=500::Int,n_sigma=2::Real)
    sample=abs.([beta_sampling(D) for i in 1:M])
    omega_m,omega_std=mean(sample),std(sample)
    T=2pi/(omega_m+n_sigma*omega_std)
    return collect(0:T/len:T)
end

"""
...
# Arguments
- `t::AbstractVector{<:Real}`: the time array for the decay curve.
- `n_D::Integer`: the number of samplings on D set.
- `sampling_D`: the function to sample over D
...
"""
function ensemble_average_FID(
    t::AbstractVector{<:Real}, 
    n_D::Int,
    sampling_D)
    
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:n_D
        f_d=ensemble_FID(t, sampling_D())
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    return f_sum/n_D, f_var/(n_D-1)
end

end