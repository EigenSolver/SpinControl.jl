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
    @time T = rabiperiod(ensemble, h, M=M, N=N, λ=0.02, L=20)
    println("period: $T")
    vals=rabi([T/2, T], ensemble,h, M=M,N=N)
    err=abs(vals[1])
    println("Error:", err)
    println("Rabi pi: ", vals[2])
    @test err<1/sqrt(M)
end