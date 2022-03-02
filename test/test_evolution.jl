

ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
cluster = SpinCluster(ensemble)

@testset "test evolution" begin
    ψ=[1,1]/√2
    ρ=ψ * ψ'
    h=10; t=π/h;
    P=SquarePulse(h, t, [1/√3, 1/√3, 1/√3])
    U=rotation(P)

    @test isunitary(U)
    @test abs(tr(ρ)-1)<1e-5
    @test evolution(ψ,U)==U*ψ
    @test evolution(ρ,U)==U*ρ*U'
    @test abs(tr(U*ρ*U')-1)<1e-5
end


@testset "test operation" begin
    ψ=[1,0]
    ρ=ψ * ψ'
    h=10; t=π/h; aim=[1,0,0]

    ϕ_p, n_p = driving(h, t, cluster, aim, N=500, sampling=true)
    P = operation(ρ, ϕ_p, n_p)
    println("trace: ", tr(P))
    @test isunitary(P)==false
    println("fidelity: ",statefidelity(P,ρ))
end
