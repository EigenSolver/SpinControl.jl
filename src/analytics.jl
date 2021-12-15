
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
function analyticalrabi(t::AbstractVector, h::Real, b::Real; axis=3)
    @assert axis in (1,2,3)
    A=(1.0.+b^4*t.^2/h^2).^(-1/4)/2
    φ=atan.(b^2*t/h)/2
    if axis==3
        return A.*cos.(h*t+φ)
    elseif axis==2
        return -A.*sin.(h*t+φ)
    end
end

