
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
# this is faster
function processfidelity(U::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Real
    F=entanglementfidelity(U, krausops)
    d=2
    return (d*F+1)/(d+1)
end

# function processfidelity(U::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Real
#     @assert isunitary(U) 
#     F=1/4
#     for ρ_0 in [σ_x, σ_y, σ_z]
#         F+=tr(U*ρ_0'*U' * operate(ρ_0,krausops))/8
#     end
#     return F
# end

function entanglementfidelity(U::Matrix{<:Number}, krausops::AbstractVector{<:Matrix})::Real
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
Return the entanglement fidelity of pauli gates (X,Y,Z) for given driving parameters
"""
# add error !!
function paulifidelity(ϕ_k::Vector{<:Real}, n_k::Matrix{<:Real}, 
    c_k::Vector{<:Real}=normalize!(ones(size(ϕ_k)),1); axis::Int = 1, phase::Symbol=:_180,
    geterr::Bool=false)
    
    if phase==:_180
        F_c = (n_k[:, axis].*sin.(ϕ_k./2)).^2 

    elseif phase==:_90
        F_c = 1/2 .* (cos.(ϕ_k./2) - n_k[:, axis].*sin.(ϕ_k./2)).^2 
    else
        error("invalid phase")
    end

    F = c_k .* F_c |>sum 
    if geterr
        F_std = sqrt.(c_k .* (F_c .- F).^2) ./ length(c_k[c_k .!=0 ])
        return F, F_std
    else
        return F
    end
end

function carrfidelity(β::AbstractVector{<:Real}, h::Real, τ::Real,
    c::Vector{<:Real}=normalize!(ones(size(β)),1); ϵ::Real=0)
    N=length(β)
    F=0
    for i in 1:N
        F+=c[i]*carrfidelity(β[i],h,τ,ϵ)
    end
    return F
end

function carrfidelity(b::Real, h::Real, τ::Real, ϵ::Real=0)::Real
    Δ=π/h+ϵ
    return 1/(4*(h^2 + b^2)^2)*(h^2 - h^2*cos(2*b*τ) 
    - cos(sqrt(h^2 + b^2)*Δ)*(h^2 + (h^2 +2*b^2)*cos(2*b*τ)) 
    + 2*b*sqrt(h^2 + b^2)*sin(sqrt(h^2 + b^2)*Δ)*sin(2*b*τ))^2
end

function xyfidelity(β::AbstractVector{<:Real}, h::Real, τ::Real,
    c::Vector{<:Real}=normalize!(ones(size(β)),1); ϵ::Real=0)
    N=length(β)
    F=0
    for i in 1:N
        F+=c[i]*xyfidelity(β[i],h,τ,ϵ)
    end
    return F
end


function xyfidelity(b::Real, h::Real, τ::Real, ϵ::Real=0)::Real
    Δ=π/h+ϵ
    Ω=sqrt(b^2+h^2)
    F=1/(64*(h^2 + b^2)^4)*
    (-5*h^4 - 8*h^2*b^2 + 3*h^4*cos(4*b*τ) + 
    8*h^2*(h^2 + 2*b^2)*cos(2*b*τ)^2 *cos(Δ*Ω) + 
    (h^4+ (h^4 + 8*h^2*b^2 + 8* b^4)*cos(4*b*τ)) *cos(2*Δ*Ω) - 
    2*h^2 *b*Ω*cos(4*b*τ - 2*Δ*Ω) - 4*b^3*Ω*cos(4*b*τ - 2*Δ*Ω) - 
    4*h^2 *b*Ω*cos(4*b*τ - Δ*Ω) + 4*h^2*b*Ω*cos(4*b*τ + Δ*Ω) + 
    2*h^2 *b*Ω*cos(4*b*τ + 2*Δ*Ω) + 4*b^3*Ω*cos(4*b*τ + 2*Δ*Ω))^2
    return F
end