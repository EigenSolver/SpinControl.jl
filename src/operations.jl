
σ_i = [1 0; 0 1]

σ_x = [0 1; 1 0]
σ_y = [0 -1im; 1im 0]
σ_z = [1 0; 0 -1]
σ_vec = (σ_x, σ_y, σ_z)

Ψ = 1/√2 .*[1,0,0,1]

import LinearAlgebra:I, tr, normalize!, normalize, norm, kron

isunitary(U::Matrix{<:Number})::Bool =  norm(U' * U- I)<1e-6

"""
Get the rotation unitary for given driving phase and axis
"""
function rotation(phi::Real, n::Vector{<:Real})::Matrix{<:Number}
    # @assert abs(norm(n)-1)<1e-4
    σ_n = mapreduce(i -> n[i] * σ_vec[i], .+, 1:3)
    return cos(phi/2)*σ_i - 1im*sin(phi/2)*σ_n # \Oemga t/2 !!!
end

"""
Get the ideal rotation unitary of a square pulse
"""
rotation(pulse::SquarePulse)=rotation(pulse.phi, pulse.aim)


"""
Evolve a quantum state or density matrix for given unitary
"""
function evolution(ψ::Vector{<:Number}, U::Matrix{<:Number})::Vector{<:Number}
    # @assert isunitary(U)
    # @assert abs(1-norm(ψ))<1e6
    return U*ψ
end

function evolution(ψ::Vector{<:Number}, phi::Real, n::Vector{<:Real})::Vector{<:Number}
    return evolution(ψ, rotation(phi, n))
end

function evolution(ρ::Matrix{<:Number}, U::Matrix{<:Number})::Matrix{<:Number}
    # @assert isunitary(U)
    # @assert abs(tr(ρ)-1)<1e6
    return U*ρ*U'
end

function evolution(ρ::Matrix{<:Number}, phi::Real, n::Vector{<:Real})::Matrix{<:Number}
    return evolution(ρ, rotation(phi, n))
end

"""
Apply quantum operation on density state for given Kraus operators
"""
function operation(ρ::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Matrix{<:Number}
    P=zeros(size(ρ))
    for E in krausops
        P+=E*ρ*E'
    end
    return P
end

"""
Apply quantum operation on density state for given rotation unitarys
"""
function operation(ρ::Matrix{<:Number}, ϕ_k::Vector{<:Real}, n_k::Matrix{<:Real}, 
    c_k::Vector{<:Real}=normalize!(ones(size(ϕ_k)), 1) 
    )::Matrix{<:Number}

    c_k=normalize(c_k,1)
    return mapreduce(i-> c_k[i]*evolution(ρ, ϕ_k[i], n_k[i, :]), .+, 1:length(c_k))
end


"""
Given a sampled list of driving strengths and phases, return a list of Kraus operators 
"""
function krausoperators(ϕ_k::Vector{<:Real}, n_k::Matrix{<:Real}, 
    c_k::Vector{<:Real}=normalize!(ones(size(ϕ_k)), 1) 
    )::Vector{<:Matrix}
    
    return [sqrt(c_k[i])*rotation(ϕ_k[i], n_k[i, :]) for i in 1:length(c_k)]
end

"""
Calculate the fidelity between two quantum states, U1 and U2
"""
function statefidelity(ρ₁::Matrix{<:Number}, ρ₂::Matrix{<:Number})
    @assert abs(tr(ρ₁)-1)<1e6 || abs(tr(ρ₂)-1)<1e6
    ρₛ = sqrt(ρ₁)
    return tr(sqrt(ρₛ*ρ₂*ρₛ))^2 
end

function statefidelity(ψ₁::Union{Matrix,Vector}, ψ₂::Union{Matrix,Vector},)::Real
    if ψ₁ isa Vector 
        ψ₁=ψ₁*ψ₁'
    end
    if ψ₂ isa Vector 
        ψ₂=ψ₂*ψ₂'
    end
    return statefidelity(ψ₂,ψ₁)
end

"""
Calculate the fidelity between two unitaries, U1 and U2
"""
function processfidelity(U₁::Matrix{<:Number},U₂::Matrix{<:Number})::Real
    @assert isunitary(U₁) || isunitary(U₂)
    return tr(U₁' * U₂)/2
end


"""
Calculate the process fidelity between the target unitary and a quantum channel defined by given kraus operators, 
"""
function processfidelity(U::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Real
    @assert isunitary(U) 
    ρ_0 = Ψ  * Ψ'

    F=0
    for E in krausops
        T = U' * E
        F += Ψ' * kron(σ_i,T) * ρ_0 * kron(σ_i,T') * Ψ 
    end
    
    return F
end

"""
Return the fidelity of pauli gates (X,Y,Z) for given driving parameters
"""
function paulifidelity(ϕ_k::Vector{<:Real}, n_k::Matrix{<:Real}, 
    c_k::Vector{<:Real}=normalize!(ones(size(ϕ_k)),1), axis::Int = 1)::Real

    return c_k .* (n_k[:, axis].*sin.(ϕ_k./2)).^2 |> sum
end

