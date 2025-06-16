ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
cluster = SpinCluster(ensemble)

# @testset "test pi fidelity" begin
#     h=50; t=π/h; 
    
#     aim=[1,0,0]
#     ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
#     U = unitary(SquarePulse(h,t,aim))

#     F1=entanglementfidelity(U, krausoperators(ϕ_p, n_p))
#     F2=paulifidelity(ϕ_p, n_p)
#     println("X gate")
#     println("Fidelity is F1:$(F1), F2:$(F2)")
#     @test abs(F1-F2)<1e-6


#     aim=[0,1,0]
#     ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
#     U = unitary(SquarePulse(h,t,aim))

#     F1=entanglementfidelity(U, krausoperators(ϕ_p, n_p))
#     F2=paulifidelity(ϕ_p, n_p, axis=2)
#     println("Y gate")
#     println("Fidelity is F1:$(F1), F2:$(F2)")
#     @test abs(F1-F2)<1e-6


#     aim=[0,0,1]
#     ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
#     U = unitary(SquarePulse(h,t,aim))

#     F1=entanglementfidelity(U, krausoperators(ϕ_p, n_p))
#     F2=paulifidelity(ϕ_p, n_p, axis=3)
#     println("Z gate")
#     println("Fidelity is F1:$(F1), F2:$(F2)")
#     @test abs(F1-F2)<1e-6
# end

@testset "test pi/2 fidelity" begin
    h=50; t=π/(2h); 
    
    aim=[-1,0,0]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = unitary(SquarePulse(h,t,aim))

    F1=entanglementfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=1, phase=:_90)
    println("√X gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6


    aim=[0,-1,0]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = unitary(SquarePulse(h,t,aim))

    F1=entanglementfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=2, phase=:_90)
    println("√Y gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6


    aim=[0,0,-1]
    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    U = unitary(SquarePulse(h,t,aim))

    F1=entanglementfidelity(U, krausoperators(ϕ_p, n_p))
    F2=paulifidelity(ϕ_p, n_p, axis=3, phase=:_90)
    println("√Z gate")
    println("Fidelity is F1:$(F1), F2:$(F2)")
    @test abs(F1-F2)<1e-6
end

@testset "test sequence fidelity" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    β=betasampling(ensemble, M=50,N=10)
    seq1=XY(20,2)
    kops=krausoperators(seq1,β)
    @time f1=processfidelity([1 0; 0 1], kops)
    @time f2=entanglementfidelity([1 0; 0 1], kops)
    @test abs(f1-(2*f2+1)/3)<1e-5

    seq2=CP(20,2)
    kops=krausoperators(seq2,β)
    @time f1=carrfidelity(β,20,2)
    @time f2=entanglementfidelity([1 0; 0 1], kops)
    println("CP fidelity:", f1)
    println("CP fidelity:", f2)
    @test abs(f1-f2)<1e-5

    seq3=XY(20,2)
    kops=krausoperators(seq3,β)
    @time f1=xyfidelity(β,20,2)
    @time f2=entanglementfidelity([1 0; 0 1], kops)
    println("XY fidelity:", f1)
    println("XY fidelity:", f2)
    @test abs(f1-f2)<1e-5
end