using Test
using SpinControl
using LinearAlgebra

# include("test_locations.jl")
# include("test_datatypes.jl")
# include("test_dynamics.jl")
# include("test_driving.jl")
# include("test_rotations.jl")
# include("test_fidelity.jl")


@testset "test period finding" begin
    M=1000; N=500; h=200;
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    T = rabiperiod(ensemble, h, M=M, N=N, λ=0.05, L=40)
    t = 0:π/h/100:π/h
    err=abs(rabi([T/2], ensemble,h, M=M,N=N)[1])
    println(err)
    @test err<1/sqrt(M)
end