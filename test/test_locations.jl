import LinearAlgebra: norm
import SpinControl.SpinEnsembles: randlocs, randsphericallocs, randcartesianlocs, randcoefs, dipolarcoefs

@testset "location generation" begin
    M=randsphericallocs(1000,1,10)

    ## Test the datatype
    @test typeof(M)<:Matrix{Float64}
    @test typeof(dipolarcoefs(M))<:Vector{Float64}
    @test all(x->1<abs(norm(x))<10, eachrow(M))

    M=randcartesianlocs(1000,1,10)
    @test typeof(M)<:Matrix{Float64}
    @test all(v->any(x->x>1,abs.(v)), eachrow(M))
    @test all(v->all(x->x<10,abs.(v)), eachrow(M))

    ## Type pass for all the function
    @test dipolarcoefs(M) isa Vector{Float64}
end 