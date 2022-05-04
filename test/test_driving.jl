ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
cluster = SpinCluster(ensemble)
h=10; t=π/h;

@testset "test pulse rotation" begin
    P=SquarePulse(h, t)
    U=rotation(P)
    @test U isa Matrix
    @test isunitary(U)

    P=SquarePulse(h, t, [1/√3, 1/√3, 1/√3])
    U=rotation(P)
    @test U isa Matrix
    println("norm: ", norm(U * U' - I))
    @test isunitary(U)
end

@testset "test rabi sampling" begin
    ω_p, n_p = rabisampling(h, t, cluster, N=500)
    ω_p2=[norm(vec) for vec in eachrow(n_p)]
    @test sum(ω_p2.-2*ω_p)<0.001
    @test [abs(norm(vec)-1)<0.001 for vec in eachrow(n_p)]|>all

    for i in 1:10
        U_p=rotation(0.2,n_p[i,:])
        @test U_p isa Matrix
        @test isunitary(U_p)
    end
end 

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