using LinearAlgebra
using QuantumOptics
using Statistics

"""
Args:
    N: number of point generated
    dim: dimension of the space
    a: scaling factor 
Return:
    an array of N random location vectors, in given dimension, scaled by a
"""
rand_loc(N::Int,dim::Int,a=1.0::Float64)::Array{Vector{Float64}}=a*collect(rand(Float64,dim) for i in 1:N)

"""
Args:
    r: point vectors from one loc to another  
    vec1,vec2:location vectors
Return:
    D: coefficient of the dipolar interaction
"""
function dipolar_coef(r::Vector{Float64})
    cosθ=r[end]/norm(r) #calculate the cos(θ) between the vector and the z axis
    D=0.5(1-3cosθ^2)/norm(r)^3
    return D
end;

dipolar_coef(vec1::Vector{Float64},vec2::Vector{Float64})=dipolar_coef(vec1-vec2);

"""
Args:
    vec0: the location of the center spin
    vec_bath: the location set of the spins in environment
Returns:
    the list of dipolar coupling strength between the centered spin and bath
"""
dipolar_bath_coefs(vec0::Vector{Float64},vec_bath::Array{Vector{Float64}})=map(x->dipolar_coef(vec0,x),vec_bath);
dipolar_bath_coefs(vec_bath::Array{Vector{Float64}})=map(dipolar_coef, vec_bath);


"""
Randomly distribute a spin bath, generate the dipolar coupling strength between the centered spin and bath
================
Args:
    [p0]: loc of center spin
    N: number of spin in bath 
    dim: dimension
    a: scale of ensemble
"""
rand_dipolar_bath_coefs(p0::Vector{Float64},N::Int,dim::Int,a=1.0::Float64)=dipolar_bath_coefs(p0,rand_loc(N,dim,a))
# place the spin at center if the p0 is not given
rand_dipolar_bath_coefs(N::Int,dim::Int,a=1.0::Float64)=dipolar_bath_coefs(collect(repeat([a/2], dim)),rand_loc(N,dim,a))

"""
Args:
    t: time, float or array of float
    D: coupling strength of dipolar interactions, usually a array of floats
Returns:
    ensemble free induction decay curve at given time t
"""
ensemble_FID(t::Float64,D::Vector{Float64})=mapreduce(cos,*,D*t)/2;
ensemble_FID(t::Vector{Float64},D::Vector{Float64})=map(x->ensemble_FID(x,D),t);


"""
Find the amplitude of tranverse magnetic field being stronger than given threshold of the coupling strength
=======================
Args:
    f: density of the spins
    a: 1d scale of the ensemble
    p: probability threshold
"""
function transverse_threshold(p::Float64, f::Float64, d::Int, a::Float64)
    N=floor(Int,a^d*f)
    D=rand_dipolar_bath_coefs(N,d,a)
    quantile(D, p)
end;