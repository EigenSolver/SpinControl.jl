@doc raw"""
    betasampling(D; N)

Give a random sampling of beta, which is a combination of `D_j`,
```math
\beta_p = \sum_j p_j \,D_j,\; p_j=\pm 1
```

# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
- `ensemble`: a `SpinEnsemble`, use this function as the method of the type 
# Options
- `N`: size of Monte-Carlo sampling
"""
function betasampling(D::Vector{<:Real}; N=1::Int)
    n=length(D)
    return [sum(rand([1,-1],n).*D) for i in 1:N]
end


@doc raw"""
    dephasingtime(D, n_t=500; scale=1.0)

The FID is Fourier transform of the noise spectrum. 
For a Gaussian noise with linewidth `b`, it's characteristic function is 
```math
f(t)=\int_{\infty}^\infty P(\beta) e^{-i \beta t} d\,\beta=\exp(-b^2 \,t^2/2), \quad P(\beta)=\frac{1}{\sqrt{2\pi b^2}} \exp(-\frac{\beta^2}{2b^2})
```
Thus the decay time is given by 
```math
T_2=\frac{\pi}{b}
```

# Arguments
- `D::Vector{Real}`: a set of coupling strengths
- `n_t`: size of the generated time array 
# Options 
- `scale`: scale factor to extend the T_2 
"""
function dephasingtime(D::Vector{<:Real}, n_t=500::Int; scale=1.0::Real)
    b=dipolarlinewidth(D)
    t=π/b*scale 
    dt=t/n_t
    return 0:dt:t 
end

dephasingtime(D::Vector{<:Real})=π/dipolarlinewidth(D)

@doc raw"""
    dipolarlinewidth(D)

Get the linewidth of D, which follows Gaussian distribution. 

```math
b=\sqrt{\sum_j D_j^2}
```
# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
"""
dipolarlinewidth(D::Vector{<:Real})=sqrt(mapreduce(abs2,+,D))

