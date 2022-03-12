@doc raw"""
    rabi(t, h, b; [axis])

Get analytical solution of Rabi oscillation under Gaussian noise with linewidth `b`

```math
M_z(t)=A \cos(h \,t+\varphi), \quad M_y(t)=-A \sin(h \,t+\varphi)
```
```math
A=\frac{1}{2}\left(1+b^4t^2/h^2\right)^{1/4},\quad \varphi=\frac{1}{2}\arctan(b^2 t/h)
```

# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `axis::Int`: , 1,2,3, representing x,y,z axis, set to 3 by default 
"""
function analyticalrabi(t::AbstractVector, cluster::SpinCluster, h::Real; axis::Int = 3)

    @assert axis in (1, 2, 3)
    b = dipolarlinewidth(cluster)

    F = exp.(-im * h * t) ./ sqrt.(im * b^2 / h * t .+ 1)

    if axis == 3
        return 1 / 2 * (b^2 / (b^2 + h^2) .+ h^2 / (b^2 + h^2) .* real(F))
    elseif axis == 2
        return 1 / 2 * h / sqrt(b^2 + h^2) .* imag(F)
    elseif axis == 1
        return 1 / 2 * b * h / (b^2 + h^2) * (-real(F) .+ 1)
    end
end

function analyticalrabi(
    t::AbstractVector{Float64},
    ensemble::SpinEnsemble,
    h::Real;
    axis::Int = 3,
)
    @assert axis in (1, 2, 3)
    @assert isdilute(ensemble)

    Γ = dipolarlinewidth(ensemble)
    F =
        exp.(im * (h .* t - Γ^2 * t ./ (2 * h))) .*
        erfc.(Γ * sqrt.(t ./ (2 * h)) * (1 - im) / sqrt(2))

    if axis == 3
        return 1 / 2 * (Γ^2 / (Γ^2 + h^2) .+ h^2 / (Γ^2 + h^2) .* real(F))
    elseif axis == 2
        return -1 / 2 * h / sqrt(Γ^2 + h^2) .* imag(F)
    elseif axis == 1
        return 1 / 2 * Γ * h / (Γ^2 + h^2) * (-real(F) .+ 1)
    end
end

function analyticalfid(
    t::AbstractVector{Float64},
    ensemble::SpinEnsemble
)
    @assert isdilute(ensemble)

    dim = ensemble.dim
    T2 = coherencetime(ensemble)
    return 1 / 2 * exp.(-(abs.(t / T2) .^ (dim / 3)))
end
