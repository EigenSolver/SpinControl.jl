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

@testset "test driving sampling" begin
    ω_p, n_p = driving(h, t, cluster, N=500, sampling=true)
    ω_p2=[norm(vec) for vec in eachrow(n_p)]
    @test sum(ω_p2.-2*ω_p)<0.001
    @test [abs(norm(vec)-1)<0.001 for vec in eachrow(n_p)]|>all

    for i in 1:10
        U_p=rotation(0.2,n_p[i,:])
        @test U_p isa Matrix
        @test isunitary(U_p)
    end

    ϕ, n = driving(h, t, cluster, N=500, sampling=false)
    println("omega: ", ϕ, " axis: ", n, "\n")
    @test n isa Vector
    @test length(n)==3

    U=rotation(ϕ, n)
    println("rotation: \n", U)
    @test U isa Matrix
    @test isunitary(U)
end 