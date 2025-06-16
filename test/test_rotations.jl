

ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
cluster = SpinCluster(ensemble)

@testset "test evolution" begin
    ψ=[1,1]/√2
    ρ=ψ * ψ'
    h=10; t=π/h;
    P=SquarePulse(h, t, [1/√3, 1/√3, 1/√3])
    U=unitary(P)

    @test isunitary(U)
    @test abs(tr(ρ)-1)<1e-5
    @test abs(tr(U*ρ*U')-1)<1e-5
end

@testset "pauli gate" begin
    X=-im*rotation(π, [-1, 0, 0])|>real
    Y=-im*rotation(π, [0, -1, 0])
    Z=-im*rotation(π, [0, 0, -1])|>real
    
    @test isunitary(X)
    println("X unitary: $X")

    @test isunitary(Y)
    println("Y unitary: $Y")

    @test isunitary(Z)
    println("Z unitary: $Z")

end


@testset "test operation" begin
    ψ=[1,0]
    ρ=ψ * ψ'
    h=10; t=π/h; aim=[1,0,0]

    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    @time P = operate(ρ, ϕ_p, n_p)
    @time P = operate(ρ, ϕ_p, n_p)

    println("trace: ", tr(P))
    @test isunitary(P)==false
    println("fidelity: ",statefidelity(P,ρ))
end
