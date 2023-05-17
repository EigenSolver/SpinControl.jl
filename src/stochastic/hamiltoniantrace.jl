
struct HamiltonianSeries
    t::AbstractVector{<:Real}
    H::AbstractMatrix{<:Real}
end


function evolve(ρ0::Matrix{<:Number}, H_t::HamiltonianSeries)
    ρ=ρ0
    for h in eachrow(H_t.H)
        U=unitary()
    end
end