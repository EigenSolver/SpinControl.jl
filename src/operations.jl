
σ_i = [1 0; 0 1]

σ_x = [0 1; 1 0]
σ_y = [0 -1im; 1im 0]
σ_z = [1 0; 0 -1]
σ_vec = (σ_x, σ_y, σ_z)

Ψ = 1/√2 .*[1,0,0,1]


## unitarys
isunitary(U::Matrix{<:Number})::Bool =  norm(U' * U- I)<1e-6


"""
Evolve a quantum state or density matrix for given unitary
"""
function evolve(ψ::Vector{<:Number}, U::AbstractMatrix{<:Number})
    return U*ψ
end

function evolve(ρ::Matrix{<:Number}, U::AbstractMatrix{<:Number})
    # @assert isunitary(U)
    # @assert abs(tr(ρ)-1)<1e6
    return U*ρ*U'
end


"""
Get the rotation unitary for given driving phase and axis
"""
function rotation(phi::Real, n::Vector{<:Real})::Matrix{ComplexF64}
    # @assert abs(norm(n)-1)<1e-4
    σ_n = mapreduce(i -> n[i] * σ_vec[i], .+, 1:3)
    return cos(phi/2)*σ_i - 1im*sin(phi/2)*σ_n # \Oemga t/2 !!!
end

"""
Get the rotation unitary for given driving vector and time
"""
function rotation(h::Vector{<:Real}, t::Real)::Matrix{ComplexF64}
    Ω=norm(h); 
    if Ω==0
        return σ_i
    else
        phi=Ω*t; n=h/Ω;
        return rotation(phi,n)
    end
end

## TPCP
"""
Apply quantum operation on density state for given Kraus operators
"""
function operate(ρ::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Matrix{ComplexF64}
    P=zeros(size(ρ))
    for E in krausops
        P+=E*ρ*E'
    end
    return P
end

function operate(ρ::Matrix{<:Number}, krausops::Vector{Adjoint{ComplexF64, Matrix{ComplexF64}}})::Matrix{ComplexF64}
    P=zeros(size(ρ))
    for E in krausops
        P+=E*ρ*E'
    end
    return P
end

"""
Given a sampled list of driving strengths and phases, return a list of Kraus operators 
"""
function krausoperators(ϕ::Vector{<:Real}, n::Matrix{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(ϕ)), 1) 
    )::Vector{<:Matrix}
    
    return [sqrt(c[i])*rotation(ϕ[i], n[i, :]) for i in eachindex(c)]
end

function krausoperators(h::Vector{<:Real},  β::AbstractVector{<:Real},
    c::Vector{<:Real}=normalize!(ones(size(β)), 1), z0::Vector=[0,0,1] )
    return krausoperators(rabisampling(h, β, z0)..., c)
end

"""
Apply quantum operation on density state for given rotaion unitarys
"""
function operate(ρ::Matrix{<:Number}, ϕ::Vector{<:Real}, n::Matrix{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(ϕ)), 1) 
    )::Matrix{ComplexF64}

    c=normalize(c,1)
    krausops=krausoperators(ϕ,n,c)
    return operate(ρ,krausops)
end

function measure(A::Matrix{<:Number}, ρ::Matrix{<:Number})::Real
    return tr(A*ρ)|>real
end

function measure(A::Matrix{<:Number}, ψ::Vector{<:Number})::Real
    return tr(ψ' *A *ψ)|>real
end