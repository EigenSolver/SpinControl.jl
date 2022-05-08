
@testset "test pulse" begin
    h=10; t=π/h;
    @test rotation(0,[0,0,0])==I

    P=SquarePulse(h, t)
    U=unitary(P)
    @test U isa Matrix
    @test isunitary(U)

    P=SquarePulse(h, t, [1/√3, 1/√3, 1/√3])
    U=unitary(P)
    @test U isa Matrix
    println("norm: ", norm(U * U' - I))
    @test isunitary(U)

    h=100; t=π/h;
    X=-1im*unitary(SquarePulse(h, t, [-1, 0, 0]))
    Y=-1im*unitary(SquarePulse(h, t, [0, -1, 0]))
    Z=-1im*unitary(SquarePulse(h, t, [0, 0, -1]))
    println(X)
    println(Z)
    @test norm(X-[0 1; 1 0])<1e-5
    @test norm(Y-[0 -1im; 1im 0])<1e-5
    @test norm(Z-[1 0; 0 -1])<1e-5
end


@testset "test idle" begin
    @test Idle(3) isa Pulse
    @test Idle(3) isa Idle
    @test unitary(Idle(3),0,[0,0,1])==I
end

@testset "test noisy pulse" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    β=betasampling(ensemble, M=50,N=10)

    h=10; t=π/h;
    P=SquarePulse(h, t, [1/√3, 1/√3, 1/√3])
    @time krausoperators(P, β)

    i=Idle(2)
    @time krausoperators(i, β)
end