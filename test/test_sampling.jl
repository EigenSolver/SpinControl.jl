@testset "test operators" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    cluster = SpinCluster(ensemble)

    h=10; t=π/h;
    U=evolution(h, t, cluster, N=300; return_unitary=true)
    println("evolution: ", U)
    @test U isa Matrix
    ω, n = evolution(h, t, cluster, N=300; return_unitary=true)
    println("omega: ", ω, " axis: ", n)
    @test n isa Vector
    @test length(n)==3
end