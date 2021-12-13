##
using Test
using SpinEnsembles
import LinearAlgebra: norm

include("../src/randloctions.jl")
include("../src/visualization.jl")

##
M=randsphericallocs(1000,1,10)

## Test the datatype
@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(randcoefs(100,3,10))<:Vector{Float64}

## Test the numerical range 
@test all(x->1<abs(norm(x))<10, eachrow(M))

M=randcartesianlocs(1000,1,10)

## Test the datatype
@test typeof(M)<:Matrix{Float64}
@test typeof(M)<:AbstractArray{Float64}
@test typeof(randcoefs(100,3))<:Vector{Float64}

## Test the numerical range 
@test all(v->any(x->x>1,abs.(v)), eachrow(M))
@test all(v->all(x->x<10,abs.(v)), eachrow(M))

## Type pass for all the function
@test dipolarcoefs(M) isa Vector{Float64}

## Test FID function
D=randcoefs(2000,3,10)
t=0:0.01:20
h=0.5
fid_curve=fid(t,D,h)
@test fid_curve isa Vector{Float64}
rabi_curve=rabi(t,D,h)
@test rabi_curve isa Vector{Float64}

## Test SpinEnsemble 
spins1=SpinEnsemble(1000,3,[1,1,0],3,5, :spherical)
@test abs(spins.rho-2.436)<1e3
spins2=SpinEnsemble(1000,2,[1,1,0],3,5, :cubic)
@test abs(spins.rho-19.984)<1e3
spin3=SpinEnsemble(1000,1,[1,1,0],3,5, :spherical)
spin4=SpinEnsemble(1000,1,[1,1,0],3,5, :cubic)
@test spin3.rho==spin4.rho

@test randlocs(spins1) isa Matrix{Float64}
randlocs!(spins2) 
@test spins2.locs isa Matrix{Float64}
