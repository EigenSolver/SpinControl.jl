
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