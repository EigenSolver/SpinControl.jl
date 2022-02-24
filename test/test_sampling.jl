@testset "test operators" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    cluster = SpinCluster(ensemble)

    h=10; t=π/h;
    ω, n = driving(h, t, cluster, N=300)
    println("omega: ", ω, " axis: ", n)
    @test n isa Vector
    @test length(n)==3
    U=evolution(ω, n)
    println("evolution: ", U)
    @test U isa Matrix
end