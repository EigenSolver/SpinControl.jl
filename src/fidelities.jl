
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
    c_k::Vector{<:Real}=normalize!(ones(size(ϕ_k)),1); axis::Int = 1)::Real

    return c_k .* (n_k[:, axis].*sin.(ϕ_k./2)).^2 |> sum
end

