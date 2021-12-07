##
using Test
using SpinEnsembles
import LinearAlgebra: norm

##
M=randlocsspherical(1,10,N=1000)

## Test the datatype
@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(randcoefs(100,3))<:Vector{Float64}

## Test the numerical range 
@test all(x->1<abs(norm(x))<10, eachrow(M))

M=randlocscubic(1,10,N=1000)

## Test the datatype
@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(randcoefs(100,3))<:Vector{Float64}

## Test the numerical range 
@test all(v->any(x->x>1,abs.(v)), eachrow(M))
@test all(v->all(x->x<10,abs.(v)), eachrow(M))

## Type pass for all the function
@test typeof(dipolarcoefs(M))<:Vector{Float64}

## 
M=randlocs(1000,10,4)
@test size(M)==(1000,10)
@test typeof(M)<:Matrix{Float64}
@test maximum(M)<4

## Test FID function
D=randcoefs(2000,3,10)
t=0:0.01:20
h=0.5
rabi(t,D,h)