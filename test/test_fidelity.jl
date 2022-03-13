ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
cluster = SpinCluster(ensemble)

@testset "test pi fidelity" begin
    h=50; t=π/h; 
    
    aim=[1,0,0]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = rotation(SquarePulse(h,t,aim))

    F1=processfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p)
    println("X gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6


    aim=[0,1,0]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = rotation(SquarePulse(h,t,aim))

    F1=processfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=2)
    println("Y gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6


    aim=[0,0,1]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = rotation(SquarePulse(h,t,aim))

    F1=processfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=3)
    println("Z gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6
end

@testset "test pi/2 fidelity" begin
    h=50; t=π/(2h); 
    
    aim=[-1,0,0]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = rotation(SquarePulse(h,t,aim))

    F1=processfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=1, phase=:_90)
    println("√X gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6


    aim=[0,-1,0]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = rotation(SquarePulse(h,t,aim))

    F1=processfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=2, phase=:_90)
    println("√Y gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6


    aim=[0,0,-1]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = rotation(SquarePulse(h,t,aim))

    F1=processfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=3, phase=:_90)
    println("√Z gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6
end