module PulseSequences

export Pulse, Idle, SquarePulse, SquareGates, Sequence, SquareSequence
export CP, CPMG, APCP, APCPMG, XY, YX, WAHUHA
export cycletime, cycleslice

abstract type Pulse end
abstract type Sequence end

struct Idle<:Pulse
    t::Real
end

struct SquarePulse<:Pulse
    h::Real
    t::Real
    aim::Vector{Real}
    phi::Float64

    function SquarePulse(h::Real, t::Real, aim::Vector{<:Real} = [1, 0, 0])
        @assert length(aim) == 3
        return new(h, t, aim, h*t) # h equals Rabi frequency at ideal limit!
    end
end

# TBD
struct GaussianPulse<:Pulse
end

#TBD
struct HermitePulse<:Pulse
end

## Sequence
SquareGates(h::Real, ϕ::Real=π, δ::Real=0)=[(aim=zeros(3); aim[n]=1; SquarePulse(h,ϕ/h+δ,aim)) for n in 1:3]

struct SquareSequence<:Sequence
    h::Real
    τ::Real
    δ::Real
    idle::Idle
    gates::Vector{SquarePulse}
    order::Vector{Int} #[0,-1,0,0,1,0]
    function SquareSequence(h::Real, τ::Real, δ::Real, order::Vector{Int})
        if maximum(order)>3
            gates=append!(SquareGates(h,π,δ), SquareGates(h,π/2,δ))
        else
            gates=SquareGates(h,π,δ)
        end
        return new(h, τ, δ, Idle(τ), gates, order)
    end

    function SquareSequence(h::Real, τ::Real, order::Vector{Int})
        return SquareSequence(h,τ,0, order)
    end
end 

struct GeneralSequence<:Sequence
    idle::Idle
    gates::Vector{<:Pulse}
    order::Vector{Int} #[0,-1,0,0,1,0]
end

function cycletime(seq::Sequence)::Real
    t0=seq.idle.t
    ts=[g.t for g in seq.gates]
    Δt=0
    for i in seq.order
        Δt+= i==0 ? t0 : ts[abs(i)] 
    end
    return Δt
end

function cycleslice(seq::Sequence, n::Int)::Vector{<:Real}
    t_now = 0.0; t_arr = [t_now];
    τ=seq.idle.t; dt=τ/n;
    t0=collect(dt:dt:τ)
    for i in seq.order
        if i == 0
            append!(t_arr, t_now .+ t0)
            t_now+=τ
        else
            t_g = seq.gates[abs(i)].t
            push!(t_arr, t_now+t_g)
            t_now+=t_g
        end
    end
    return t_arr
end


include("protocols.jl")

end