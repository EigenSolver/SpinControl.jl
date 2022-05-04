module PulseSequences

export Pulse, Idle, SquarePulse, SquareGates
export Sequence, SquareSequence, CP, CPMG, XY, YX

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
    idle::Idle
    gates::Vector{SquarePulse}
    order::Vector{Int} #[0,-1,0,0,1,0]
    function SquareSequence(h::Real, τ::Real, δ::Real, order::Vector{Int})
        gates=SquareGates(h,π,δ)
        return new(Idle(τ), gates, order)                
    end
end 

struct GeneralSequence<:Sequence
    idle::Idle
    gates::Vector{<:Pulse}
    order::Vector{Int} #[0,-1,0,0,1,0]
end

CP(h::Real,τ::Real, δ::Real=0)=SquareSequence(h, τ, δ,[0,1,0,0,1,0])
APCP(h::Real,τ::Real, δ::Real=0)=SquareSequence(h, τ, δ, [0,-1,0,0,1,0])
CPMG(h::Real,τ::Real, δ::Real=0)=SquareSequence(h, τ, δ, [0,2,0,0,2,0])
APCPMG(h::Real,τ::Real, δ::Real=0)=SquareSequence(h, τ, δ, [0,-2,0,0,2,0])

order_xy=[0,1,0,0,2,0,0,1,0,0,2,0]

function XY(h::Real,τ::Real, δ::Real=0; symmetry::Bool=false)
    if symmetry
        return SquareSequence(h, τ, δ, vcat(order_xy, reverse!(order_xy)))
    else
        return SquareSequence(h, τ, δ, order_xy)
    end
end

function YX(h::Real,τ::Real, δ::Real=0; symmetry::Bool=false)
    if symmetry
        return SquareSequence(h, τ, δ,vcat(reverse!(order_xy),order_xy))
    else
        return SquareSequence(h, τ, δ, reverse!(order_xy))
    end
end

end