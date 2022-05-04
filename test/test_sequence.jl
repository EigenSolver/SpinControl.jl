
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


@testset "test sequence" begin
    seq1=XY(10,1)
    @test seq1 isa Sequence
    @test seq1.idle isa Idle
    @test seq1.idle isa Pulse
    @test seq1.gates isa Vector{<:Pulse}
    @test seq1.order isa Vector

    seq2=XY(10, 1, symmetry=true)
    @test seq2 isa Sequence
    println("XY-8: ", seq2)
    seq3=APCP(10,1)
    @test seq3 isa Sequence
    println("APCP: ", seq3)
    
    I1=unitary(seq1)
    I2=unitary(seq2)
    I3=unitary(seq3)
    println(I1)
    @test norm(I1+I)<1e-5
    @test norm(I2-I)<1e-5
    @test norm(I3-I)<1e-5
end

@testset "test noisy sequence" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    β=betasampling(ensemble, M=50,N=10)
    seq1=XY(100,2)
    kops=krausoperators(seq1,β)
    ρ=[1 -im; im 1]/2
    ρt=operation(ρ, kops)
    println(ρt)
    f=statefidelity(ρ,ρt)
    println(f)
end