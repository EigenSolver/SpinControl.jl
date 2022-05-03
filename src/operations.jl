
σ_i = [1 0; 0 1]

σ_x = [0 1; 1 0]
σ_y = [0 -1im; 1im 0]
σ_z = [1 0; 0 -1]
σ_vec = (σ_x, σ_y, σ_z)

Ψ = 1/√2 .*[1,0,0,1]

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
function operation(ρ::Matrix{<:Number}, ϕ::Vector{<:Real}, n::Matrix{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(ϕ)), 1) 
    )::Matrix{<:Number}

    c=normalize(c,1)
    return mapreduce(i-> c[i]*evolution(ρ, ϕ[i], n[i, :]), .+, 1:length(c))
end


"""
Given a sampled list of driving strengths and phases, return a list of Kraus operators 
"""
function krausoperators(ϕ::Vector{<:Real}, n::Matrix{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(ϕ)), 1) 
    )::Vector{<:Matrix}
    
    return [sqrt(c[i])*rotation(ϕ[i], n[i, :]) for i in 1:length(c)]
end