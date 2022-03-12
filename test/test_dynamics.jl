import Statistics: mean, std

function testcurve(
    the_curve::Vector{Float64},
    mc_curve::Vector{Float64},
    avg_bound::Real,
    std_bound::Real,
)

    error = abs.(mc_curve - the_curve)
    avg_err = mean(error)
    std_err = std(error)

    println("mean error: ", avg_err)
    println("std error: ", std_err)
    @test avg_err < avg_bound
    @test std_err < std_bound
end

@testset "cluster dynamics" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    cluster = SpinCluster(ensemble)
    N = 500
    h = 40
    T2 = coherencetime(cluster)
    t = 0:π/(20*h):T2

    avg_bound = std_bound = 2 / √N # important

    println("testing cluster rabi oscillation curve Sz...")
    testcurve(
        analyticalrabi(t, cluster, h),
        rabi(t, cluster, h; N = N),
        avg_bound,
        std_bound,
    )

    println("testing cluster rabi oscillation curve Sy...")
    testcurve(
        analyticalrabi(t, cluster, h, axis = 2),
        rabi(t, cluster, h; N = N, axis = 2),
        avg_bound,
        std_bound,
    )

    println("testing cluster rabi oscillation curve Sx...")
    testcurve(
        analyticalrabi(t, cluster, h, axis = 1),
        rabi(t, cluster, h; N = N, axis = 1),
        avg_bound,
        std_bound,
    )
end

@testset "ensemble average dynamics" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    T2 = coherencetime(ensemble)
    @test abs(T2 - 1) < 0.001
    h = 50
    dt = π / (h * 20)
    t = 0:dt:T2
    M = 600
    N = 500

    avg_bound = std_bound = 10 / sqrt(M * N)

    println("testing ensemble rabi oscillation curve Sz...")
    testcurve(
        analyticalrabi(t, ensemble, h),
        rabi(t, ensemble, h; M = M, N = N),
        avg_bound,
        std_bound,
    )


    println("testing ensemble rabi oscillation curve Sy...")
    testcurve(
        analyticalrabi(t, ensemble, h, axis = 2),
        rabi(t, ensemble, h; axis = 2, M = M, N = N),
        avg_bound,
        std_bound,
    )


    println("testing ensemble rabi oscillation curve Sx...")
    testcurve(
        analyticalrabi(t, ensemble, h, axis = 1),
        rabi(t, ensemble, h; axis = 1, M = M, N = N),
        avg_bound,
        std_bound,
    )

    dt = T2 / 200
    t = 0:dt:T2
    println("testing ensemble free induction decay curve...")
    testcurve(
        analyticalfid(t, ensemble),
        fid(t, ensemble; M = M, N = N),
        avg_bound,
        std_bound,
    )
end

@testset "test period finding" begin
    M=1000; N=500; h=200;
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    T = rabiperiod(ensemble, h, M=M,N=N)
    t = 0:π/h/100:π/h
    @test abs(rabi([T/2], ensemble,h, M=M,N=N)[1])<1/sqrt(M)
end