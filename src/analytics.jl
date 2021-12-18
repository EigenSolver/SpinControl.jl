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
function analyticalrabi(t::AbstractVector, cluster::SpinCluster, h::Real; 
    axis=3::Int, spin=1/2::Real)

    @assert axis in (1,2,3)
    b=dipolarlinewidth(cluster)
    A=(1.0.+b^4*t.^2/h^2).^(-1/4)
    φ=atan.(b^2*t/h)/2
    F=spin*A.*exp.(-im*h*t+φ)

    if axis==3
        return real(F) 
    elseif axis==2
        return imag(F)
    elseif axis==1
        return zeros(length(t))
    end
end

function analyticalrabi(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h::Real;
    axis=3::Int, spin=1/2::Real)
    @assert axis in (1,2,3)
    @assert isdilute(ensemble)

    Γ=coherencetime(ensemble)
    F=spin.* exp.(im*(h.*t-Γ^2*t ./(2*h))) .* erfc.(Γ*sqrt.(t ./(2*h)) *(1-im)/sqrt(2))

    if axis==3
        return real(F)
    elseif axis==2
        return imag(F)
    else
        return zeros(length(t))
    end
end

function analyticalfid(t::AbstractVector{Float64}, ensemble::SpinEnsemble;
     spin=1/2::Real)
    @assert isdilute(ensemble)

    dim=ensemble.dim
    T2=coherencetime(ensemble)
    return spin*exp.(-(abs.(t/T2).^(dim/3)))
end