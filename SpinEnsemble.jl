using LinearAlgebra
using Statistics
using Random

include("RandLoctions.jl")
"""
Given the locations of central spin and its bath spin, return the vector set from the central spin to bath 
===============
Args:
    loc0: location of the central spin
    loc_bath: location set of the spin bath
Return:
    the vector set from the central spin to bath 
"""
bath_vectors(loc0::Vector{<:Real},loc_bath::Matrix{<:Real})=map(x->x-loc0,eachrow(loc_bath))

"""
Calculate the dipolar interaction strength given r and z0
===============
Args:
    r: point vectors from one loc to another  
    vec1,vec2:location vectors
    z0: direction of the background magnetic field, default to be [0,0,1]
Return:
    D: coefficient of the dipolar interaction
"""
function dipolar_coef(r::AbstractArray{<:Real},z0::AbstractArray{<:Real})
    # suppose z0 is already normalized
    cosθ=dot(r,z0)/norm(r) #calculate the cos(θ) between the vector and the z axis
    D=0.5*(1-3cosθ^2)/norm(r)^3
    return D
end;


"""
Args:
    vec_bath: an array of vector, distance from the central spin to the spins in bath 
    z0: z axis
Returns:
    the list of dipolar coupling strength between the centered spin and bath
"""
function bath_dipolar_coefs(vec_bath::Matrix{<:Real},z0=[0,0,1.0]::Vector{<:Real})
    normalize!(z0)
    map(x->dipolar_coef(x,z0), eachrow(vec_bath))
end
# use for the case when the central spin is not at zero
"""
Randomly distribute a spin bath, generate the dipolar coupling strength between the centered spin and bath
===============
Args:
    N: number of spin in bath 
    dim: dimension
    a: scale of ensemble
"""
rand_bath_dipolar_coefs(N::Int,dim::Int,scaling_a=1::Real)=bath_dipolar_coefs(rand_locs(N,dim,a))

function rand_bath_dipolar_coefs(N::Int, bound::Tuple{Real,Real}, method=:shperical)
    
end

"""
Args:
    t: time, float or array of float
    D: coupling strength of dipolar interactions, usually a array of floats
Returns:
    ensemble free induction decay curve at given time t
"""
ensemble_FID(t::Real,D::Vector{<:Real})=mapreduce(cos,*,D*t)/2;
ensemble_FID(t::AbstractVector{<:Real},D::Vector{<:Real})=map(x->ensemble_FID(x,D),t);


"""
Find the amplitude of tranverse magnetic field being stronger than given threshold of the coupling strength
===============
Args:
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
===============
Args:
    D_set: a set of the coupling strengths
"""
beta_sampling(D_set::Vector{<:Real})=sum(rand([1,-1],length(D_set)).*D_set) 

"""
Get a random sampling of the f=Sx(t), under given transverse magnetic field
===============
Args:
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
===============
Args:
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
Given a set of coupling strength, determine the maximum time scale and minimum time scale required for the problem
===============
"""
function ensemble_average_FID(M::Int, N::Int)
    return 0
end


