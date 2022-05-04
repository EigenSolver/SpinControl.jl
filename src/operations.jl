
σ_i = [1 0; 0 1]

σ_x = [0 1; 1 0]
σ_y = [0 -1im; 1im 0]
σ_z = [1 0; 0 -1]
σ_vec = (σ_x, σ_y, σ_z)

Ψ = 1/√2 .*[1,0,0,1]


## unitarys
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
Get the rotation unitary for given driving vector and time
"""
function rotation(h::Vector{<:Real}, t::Real)::Matrix{<:Number}
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
function operation(ρ::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Matrix{<:Number}
    P=zeros(size(ρ))
    for E in krausops
        P+=E*ρ*E'
    end
    return P
end

"""
Apply quantum operation on density state for given rotaion unitarys
"""
function operation(ρ::Matrix{<:Number}, ϕ::Vector{<:Real}, n::Matrix{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(ϕ)), 1) 
    )::Matrix{<:Number}

    c=normalize(c,1)
    krausops=krausoperators(ϕ,n,c)
    return operation(ρ,krausops)
end

# rabisampling->krausoperators

"""
Given a sampled list of driving strengths and phases, return a list of Kraus operators 
"""
function krausoperators(ϕ::Vector{<:Real}, n::Matrix{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(ϕ)), 1) 
    )::Vector{<:Matrix}
    
    return [sqrt(c[i])*rotation(ϕ[i], n[i, :]) for i in 1:length(c)]
end

# rabisampling->paulifidelity
"""
Get the unitary of a square pulse
"""

function unitary(pulse::SquarePulse, β::Real=0, z0::Vector{<:Real}=[0,0,1])::Matrix{<:Number}
    return rotation(pulse.aim.*pulse.h + z0.*β, pulse.t)
end

function unitary(pulse::Idle, β::Real=0, z0::Vector{<:Real}=[0,0,1])
    return rotation(z0.*β, pulse.t)
end

function unitary(seq::Sequence, β::Real=0, z0::Vector{<:Real}=[0,0,1])
    U0=unitary(seq.idle,β,z0)
    Un=[unitary(g,β,z0) for g in seq.gates]
    V=σ_i
    for i in seq.order
        U= i==0 ? U0 : (sign(i)>0 ? Un[abs(i)] : Un[abs(i)]')
        V=U*V
    end
    return V
end

function krausoperators(pulse::Pulse, β::AbstractVector{<:Real},
    c::Vector{<:Real}=normalize!(ones(size(β)), 1),
    z0::Vector{<:Real}=[0,0,1]
    )::Vector{<:Matrix}
    return sqrt.(c).*[unitary(pulse, β_k, z0) for β_k in β]
end

function krausoperators(seq::Sequence, β::AbstractVector{<:Real},
    c::Vector{<:Real}=normalize!(ones(size(β)), 1),
    z0::Vector{<:Real}=[0,0,1]
    )::Vector{<:Matrix}
    
    return sqrt.(c).*[unitary(seq, β_k, z0) for β_k in β]
end
