##
using Test
import LinearAlgebra: norm
include("RandLoctions.jl")
include("SpinEnsemble.jl")

## Test the numerical range 
@test minimum(x->abs(norm(x)), eachrow(rand_locs_spherical(1,10,N=100)))>1
@test maximum(x->abs(norm(x)), eachrow(rand_locs_spherical(1,10,N=100)))<10

## Test the datatype

M=rand_locs_spherical(1,10,N=1000)

@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(rand_bath_dipolar_coefs(100,3))<:Vector{Float64}

## Type pass for all the function
@test typeof(bath_dipolar_coefs(M))<:Vector{Float64}