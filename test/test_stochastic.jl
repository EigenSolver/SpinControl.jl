include("../src/stochastic/stochasticprocess.jl")

@testset "OU Process" begin
    params=(0.007,0.01,0)
    α, β, γ=params
    σ = sqrt(β^2/(2*α))
    μ = γ
    n = 50;

    x=range(-4*σ,4*σ,100);
    ou=OrnsteinUhlenbeckNoise(n, 1, params)

    @test length(ou.u) == n+1
    @test ou isa NoiseProcess

    unitary(ou)
    dephase_curve=fid(ou)
    @test length(dephase_curve)==n+1
    
    rabi_curve=rabi(ou, 100)
    @test length(rabi_curve)==n+1
end