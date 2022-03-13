

ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
cluster = SpinCluster(ensemble)

@testset "test evolution" begin
    ψ=[1,1]/√2
    ρ=ψ * ψ'
    h=10; t=π/h;
    P=SquarePulse(h, t, [1/√3, 1/√3, 1/√3])
    U=rotation(P)

    @test isunitary(U)
    @test abs(tr(ρ)-1)<1e-5
    @test evolution(ψ,U)==U*ψ
    @test evolution(ρ,U)==U*ρ*U'
    @test abs(tr(U*ρ*U')-1)<1e-5
end

@testset "pauli gate via evolution" begin
    h=100; t=π/h;
    X=rotation(SquarePulse(h, t, [-1, 0, 0]))
    Y=rotation(SquarePulse(h, t, [0, -1, 0]))
    Z=rotation(SquarePulse(h, t, [0, 0, -1]))

    X=-im*rotation(π, [-1, 0, 0])|>real
    Y=-im*rotation(π, [0, -1, 0])
    Z=-im*rotation(π, [0, 0, -1])|>real
    
    @test isunitary(X)
    println("X rotation: $X")

    @test isunitary(Y)
    println("Y rotation: $Y")

    @test isunitary(Z)
    println("Z rotation: $Z")

end


@testset "test operation" begin
    ψ=[1,0]
    ρ=ψ * ψ'
    h=10; t=π/h; aim=[1,0,0]

    ϕ_p, n_p = rabisampling(h, t, cluster, aim, N=1000)
    @time P = operation(ρ, ϕ_p, n_p)
    @time P = operation(ρ, ϕ_p, n_p)


    # function operation2(ρ::Matrix{<:Number}, ϕ_k::Vector{<:Real}, n_k::Matrix{<:Real}, 
    #     c_k::Vector{<:Real}=ones(size(ϕ_k)))::Matrix{<:Number}
    
    #     ops=krausoperators(ϕ_k, n_k,c_k)
    #     return operation(ρ,ops)
    # end

    # @time P2 = operation2(ρ, ϕ_p, n_p)
    # @time P2 = operation2(ρ, ϕ_p, n_p)

    # @test norm(P-P2)<1e5
    # println("dif: ", P-P2)

    println("trace: ", tr(P))
    @test isunitary(P)==false
    println("fidelity: ",statefidelity(P,ρ))
end
