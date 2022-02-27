
σ_i = [1 0; 0 1]
σ_x = [0 1; 1 0]
σ_y = [0 -1im; 1im 0]
σ_z = [1 0; 0 -1]
σ_vec = (σ_x, σ_y, σ_z)

import LinearAlgebra:I, tr, normalize!, norm

isunitary(U::Matrix{<:Number})::Bool =  norm(U' * U- I)<1e-6

"""
Calculate the fidelity between two unitaries, U1 and U2
"""
function fidelity(U1::Matrix{<:Number},U2::Matrix{<:Number})::Real
    @assert isunitary(U1) || isunitary(U2)
    return tr(U1' * U2)/2
end

"""
Get the evolution unitary for given driving phase and axis
"""
function evolution(phi::Real, n::Vector{<:Real})::Matrix
    @assert abs(norm(n)-1)<1e-4
    σ_n = mapreduce(i -> n[i] * σ_vec[i], .+, 1:3)
    return cos(phi)*σ_i - 1im*sin(phi)*σ_n
end

"""
Get the ideal evolution operator of a square pulse
"""
evolution(pulse::SquarePulse)=evolution(pulse.phi, pulse.aim)

"""
Get the noisy quantum channel by averageing the many-body unitary projection 
on the central qubit subspace, which equals a partial trace. 
"""
function channel(h::Real, t::Real, cluster::SpinCluster, 
    aim::Vector{<:Real}=[1,0,0]; N::Int=100)::Matrix
    ϕ_p, n_p = driving(h, t, cluster, aim, N=N, sampling=true)
    ξ = mapreduce(i->evolution(ϕ_p[i], n_p[i,:]), .+, 1:N) ./ N
    return ξ
end
