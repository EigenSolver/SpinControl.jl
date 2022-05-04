using Test
using SpinControl
using LinearAlgebra

# include("test_locations.jl")
# include("test_datatypes.jl")
# include("test_dynamics.jl")
# include("test_driving.jl")
# include("test_rotations.jl")
# include("test_fidelity.jl")


@testset "sequence" begin
    @test Idle(3) isa Pulse
    @test Idle(3) isa Idle
    seq1=XY(10,1)
    @test seq1 isa Sequence
    seq2=vcat(XY(10,1), YX(10,1))
    @test seq2 isa Sequence
    println("XY-8: ", seq2)
    seq3=CP(10,1)
    @test seq3 isa Sequence
    println("CP: ", seq3)
    # seq4=CPMG(10,1)
    # @test seq4 isa Sequence
    # println("CPMG: ", seq4)
end
