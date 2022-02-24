import SpinControl.SpinEnsembles: reroll!

@testset "spin ensemble" begin
    ensemble = SpinEnsemble(1000, 3, [1, 1, 0], 3, 5, :spherical)
    @test abs(ensemble.rho - 2.436) < 1e3
    @test volume(ensemble) == 4π / 3 * (5^3 - 3^3)
    @test isdilute(ensemble) == false
    @test randlocs(ensemble) isa Matrix{Float64}
    @test randcoefs(ensemble) isa Vector{Float64}

    ensemble = SpinEnsemble(19.984, 2, [1, 1, 0], 3, 5, :cubic)
    @test abs(ensemble.n - 1000) < 1e3
    @test volume(ensemble) == 4 * (5^2 - 3^2)

    ensemble3 = SpinEnsemble(1000, 1, [1, 1, 0], 0, 5, :spherical)
    ensemble4 = SpinEnsemble(2000, 1, [1, 1, 0], 0, 5, :cubic)
    @test ensemble3.rho == ensemble4.rho / 2
    @test volume(ensemble3) == volume(ensemble4) == 10
end

@testset "spin cluster" begin
    ensemble = SpinEnsemble(1000, 3, [1, 1, 0], 3, 5, :spherical)
    cluster = SpinCluster(ensemble)
    @test dipolarlinewidth(cluster) isa Real
    @test betasampling(cluster, 100) isa Vector{Float64}
    @test coherencetime(cluster) == π / dipolarlinewidth(cluster)
    reroll!(cluster)
end
