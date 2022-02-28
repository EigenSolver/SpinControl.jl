module PulseSequences

export Pulse, SquarePulse

abstract type Pulse end

struct SquarePulse<:Pulse
    h::Real
    t::Real
    aim::Vector{Real}
    phi::Float64

    function SquarePulse(h::Real, t::Real, aim::Vector{<:Real} = [1, 0, 0])
        @assert length(aim) == 3
        return new(h, t, aim, h*t)
    end
end

end
