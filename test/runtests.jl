##
using Test
import LinearAlgebra: norm

include("../src/RandLoctions.jl")
include("../src/Visualization.jl")
using SpinEnsembles: rand_bath_dipolar_coefs, bath_dipolar_coefs

##
M=rand_locs_spherical(1,10,N=1000)

## Test the datatype
@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(rand_bath_dipolar_coefs(100,3))<:Vector{Float64}

## Test the numerical range 
@test all(x->1<abs(norm(x))<10, eachrow(M))

M=rand_locs_cubic(1,10,N=1000)

## Test the datatype
@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(rand_bath_dipolar_coefs(100,3))<:Vector{Float64}

## Test the numerical range 
@test all(v->any(x->x>1,abs.(v)), eachrow(M))
@test all(v->all(x->x<10,abs.(v)), eachrow(M))

## Type pass for all the function
@test typeof(bath_dipolar_coefs(M))<:Vector{Float64}

## 
M=rand_locs(1000,10,4)
@test size(M)==(1000,10)
@test typeof(M)<:Matrix{Float64}
@test maximum(M)<4
