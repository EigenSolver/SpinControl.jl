##
using Test
using SpinEnsembles
import LinearAlgebra: norm


include("../src/visualization.jl")
## Test FID function
D=randcoefs(2000,3,10)
t=0:0.01:20
h=0.5
fid_curve=fid(t,D,h)
@test fid_curve isa Vector{Float64}
rabi_curve=rabi(t,D,h)
@test rabi_curve isa Vector{Float64}

## Test SpinEnsemble 
ensemble1=SpinEnsemble(1000,3,[1,1,0],3,5, :spherical)
@test abs(ensemble1.rho-2.436)<1e3
ensemble2=SpinEnsemble(1000,2,[1,1,0],3,5, :cubic)
@test abs(ensemble2.rho-19.984)<1e3
ensemble3=SpinEnsemble(1000,1,[1,1,0],3,5, :spherical)
ensemble4=SpinEnsemble(1000,1,[1,1,0],3,5, :cubic)
@test ensemble3.rho==ensemble4.rho

@test randlocs(ensemble1) isa Matrix{Float64}
randlocs!(ensemble2) 
@test ensemble2.locs isa Matrix{Float64}

## Test SpinCluster
spins1=SpinCluster(ensemble1)

