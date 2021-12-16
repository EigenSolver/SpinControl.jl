import LinearAlgebra: norm
include("../src/randloctions.jl")

@testset "location generation" begin
    M=randsphericallocs(1000,1,10)

    ## Test the datatype
    @test typeof(M)<:Matrix{Float64}
    @test typeof(M)<:AbstractVector{Float64}
    @test typeof(randcoefs(100,3,10))<:Vector{Float64}

    ## Test the numerical range 
    @test all(x->1<abs(norm(x))<10, eachrow(M))

    M=randcartesianlocs(1000,1,10)

    ## Test the datatype
    @test typeof(M)<:Matrix{Float64}
    @test typeof(M)<:AbstractVector{Float64}
    @test typeof(randcoefs(100,3))<:Vector{Float64}

    ## Test the numerical range 
    @test all(v->any(x->x>1,abs.(v)), eachrow(M))
    @test all(v->all(x->x<10,abs.(v)), eachrow(M))

    ## Type pass for all the function
    @test dipolarcoefs(M) isa Vector{Float64}
end 