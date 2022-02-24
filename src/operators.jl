
σ_i = [1 0; 0 1]
σ_x = [0 1; 1 0]
σ_y = [0 -1im; 1im 0]
σ_z = [1 0; 0 -1]
σ_vec = (σ_x, σ_y, σ_z)

import LinearAlgebra:I, tr, normalize!
isunitary(U::Matrix{<:Number}) = U'.U == I

function fidelity(U1::Matrix{<:Number},U2::Matrix{<:Number})
    @assert isunitary(U1) 
    @assert isunitary(U2)

    return tr(U1*U2)
end

function evolution(ϕ::Real, n::Vector{<:Real})
    normalize!(n)
    return cos.(ϕ)*σ_i+sin(ϕ)*sum([n[i]*σ_vec[i] for i in 1:3])
end

function evolution(pulse::Pulse)
    phi=pulse.phi; n=pulse.aim;
    σ_n = mapreduce(i -> n[i] * σ_vec[i], +, 1:3)
    U = cos(phi) * σ_i + sin(phi) * σ_n
    return U
end