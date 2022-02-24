module PulseSequences

struct SquarePulse
    h::Real
    t::Real
    aim::Vector{Real}
    phi::Float64
    U::Matrix{Complex}

    function SquarePulse(h::Real, t::Real, aim::Vector{<:Real} = [1, 0, 0])
        @assert length(aim) == 3
        phi = h * t / 2 # pauli Matrix
        σ_n = mapreduce(i -> aim[i] * σ_vec[i], +, 1:3)
        println(σ_n)
        println(typeof(σ_n))
        U = cos(phi) * σ_i + sin(phi) * σ_n
        return new(h, t, aim, phi, U)
    end
end

end
