# @testset "test sequence" begin
#     h=10; τ=1;
#     seq1=XY(h,τ)
#     @test seq1 isa Sequence
#     @test seq1.idle isa Idle
#     @test seq1.idle isa Pulse
#     @test seq1.gates isa Vector{<:Pulse}
#     @test seq1.order isa Vector

#     seq2=XY(h,τ, symmetry=true)
#     @test seq2 isa Sequence
#     println("XY-8: ", seq2)
#     seq3=APCP(h,τ)
#     @test seq3 isa Sequence
#     println("APCP: ", seq3)
    
#     @test cycletime(seq1)==8*τ+4π/h
#     @test cycletime(seq2)==16*τ+8π/h
#     @test cycletime(seq3)==4*τ+2π/h
    
#     @test cycleslice(seq3,2)==
#     [0,τ/2,τ,τ+π/h,1.5τ+π/h,2τ+π/h,2.5τ+π/h,3τ+π/h,3τ+2π/h,3.5τ+2π/h,4τ+2π/h]

#     I1=unitary(seq1)
#     I2=unitary(seq2)
#     I3=unitary(seq3)
#     println(I1)
#     @test norm(I1+I)<1e-5
#     @test norm(I2-I)<1e-5
#     @test norm(I3-I)<1e-5
# end

# @testset "test XY sequence" begin
#     ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
#     β=betasampling(ensemble, M=500,N=100)
#     seq1=XY(30,0.2)
#     println("kraus ops time:")
#     @time kops=krausoperators(seq1,β)
#     ρ=[1 -im; im 1]/2

#     println("operation time:")
#     ρt=operate(ρ, kops)
#     # println(ρt)
#     f=statefidelity(ρ,ρt)
#     println("fidelity:", f)


#     println("trace time:")
#     ψ=[1, im]/sqrt(2)
#     @time t_arr, ψ_arr= deploy(ψ, seq1, 5, β[1])
#     @test length(t_arr)==length(ψ_arr)
#     @time t_arr, ψ_arr= deploy(ψ, seq1, 5, β[1], cycle=20)
#     @test length(t_arr)*20==length(ψ_arr)

#     @time t_arr, ρ_arr= deploy(ρ, seq1, 5, β)
#     @test length(t_arr)==length(ρ_arr)
#     @time t_arr, ρ_arr= deploy(ρ, seq1, 5, β, cycle=2)
#     @test length(t_arr)*2==length(ρ_arr)
    
# end

@testset "test WAHUHA sequence" begin
    ensemble = SpinEnsemble(0.39486, 3, [0, 0, 1], 0.1, 10, :spherical)
    β=betasampling(ensemble, M=500,N=100)
    seq1=WAHUHA(30,0.2)
    println("kraus ops time:")
    @time kops=krausoperators(seq1,β)
    ρ=[1 -im; im 1]/2

    println("operation time:")
    ρt=operate(ρ, kops)
    # println(ρt)
    f=statefidelity(ρ,ρt)
    println("fidelity:", f)

    println("trace time:")
    ψ=[1, im]/sqrt(2)
    @time t_arr, ψ_arr= deploy(ψ, seq1, 5, β[1])
    @test length(t_arr)==length(ψ_arr)
    @time t_arr, ψ_arr= deploy(ψ, seq1, 5, β[1], cycle=20)
    @test length(t_arr)*20==length(ψ_arr)

    @time t_arr, ρ_arr= deploy(ρ, seq1, 5, β)
    @test length(t_arr)==length(ρ_arr)
    @time t_arr, ρ_arr= deploy(ρ, seq1, 5, β, cycle=2)
    @test length(t_arr)*2==length(ρ_arr)
    
end