

@doc raw"""
    fid(t, D)

Given a set of dipolar coupling constant, 
calculate the free induction decay (FID) value at given time `t`.

```math
f(t)=\frac{1}{2}\prod_j \cos(D_jt)
```
# Arguments
- `t`: time, float or array of float
- `D`: coupling strength of dipolar interactions, a array of floats
"""
fid(t::Real,D::Vector{<:Real})=mapreduce(cos,*,D*t)/2;
fid(t::AbstractVector{<:Real},D::Vector{<:Real})=map(x->fid(x,D),t);

"""
    fid(spins; options)

Calculate the free induction dacay of the given spin cluster
"""
function fid(spins::SpinCluster; n_t=200::Int, scale=1.2::Real)
    t=dephasingtime(spins)*scale
    return fid(0:t/n_t:t, spins.couplings)
end

@doc raw"""
    averagefid(t, M, sampling_D; [options]...)
    averagefid(t, M, sampling_D; [options]...)

Calculate the average free induction decay over different ensembles (disorders) 

```math
\bar{f}(t)=\sum_k f(t; \{D_i\}_k)
```

# Arguments
- `t`: the time array for the decay curve.
- `M`: number of ensembles  
- `ensemble`: spin ensemble
- `sampling_D`: the function to sample over D
"""
function averagefid(t::AbstractVector{<:Real}, M::Int, sampling_D::Function; returnerr=true)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=fid(t, sampling_D())
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
end

function averagefid(ensemble::SpinEnsemble; M=1000::Int, n_t=200::Int, scale=1.2::Real)
    t=dephasingtime(ensemble, n_t; scale=scale)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=fid(t, randcoefs(ensemble))
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
end


@doc raw"""
    fid(t, D, h; N)
    fid(t, spins, h; N)

Simulate FID signal with given transverse magnetic field `h`, using Monte-Carlo sampling

```math
f_p(t)=\frac{1}{2}[\cos^2(\omega_p t)+\sin^2(\omega_p t) (n_x^2-n_z^2)]
```

```math
\bar{f}(t)=\frac{1}{N}\sum_{p=1}^N f_p(t)
```
# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `N`: number of Monte-Carlo sampling
- `geterr::Bool`: wetehr to return the error of the monte-sampling 
"""
function fid(t::AbstractVector{<:Real},D::Vector{<:Real},h::Real; N=100::Int, geterr=:false)
    n=length(D)
    f_sum=zeros(length(t)) # sum
    f_var=copy(f_sum) # square sum
    for i in 1:N
        β_p=sum(rand([1,-1],n).*D) 
        ω_p=sqrt(h^2+β_p^2)/2
        cos2_p=cos.(ω_p*t).^2
        f_p=(cos2_p+(cos2_p.-1)*(β_p^2-h^2)/(h^2+β_p^2))/2
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum).^2/(i*(i-1)) : f_var
    end
    if geterr
        return (f_sum/N,f_var/(N-1))
    else
        return f_sum/N
    end
end

fid(t::AbstractVector{<:Real},spins::SpinCluster,h::Real; N=100::Int, geterr=:false)=
fid(t, spins.couplings, h; N=N, geterr=geterr)

function fid(spins::SpinCluster,h::Real; 
    N=100::Int, n_t=200::Int, scale=1.2::Real, geterr=:false)
    t=dephasingtime(spins)*scale
    fid(0:t/n_t:t, spins.couplings, h; N=N, geterr=geterr)
end


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
function rabi(t::AbstractArray, h::Real, b::Real; axis=3)
    @assert axis in (1,2,3)
    A=(1.0.+b^4*t.^2/h^2).^(-1/4)/2
    φ=atan.(b^2*t/h)/2
    if axis==3
        return A.*cos.(h*t+φ)
    elseif axis==2
        return -A.*sin.(h*t+φ)
    end
end


@doc raw"""
    rabi(t, D, h; N, options...)

Get a random sampling of, under given transverse magnetic field

```math
\begin{aligned}
z_p(t)&=+\frac{1}{2} \left[n_z^2+n_x^2 \cos (\Omega_p t)\right]\\
y_p(t)&=-\frac{1}{2} n_x \sin (\Omega_p t)\\
x_p(t)&=+\frac{1}{2} \left[n_x n_z-n_x n_z \cos (\Omega_p t)\right]
\end{aligned}
```

```math
G(t)=\frac{1}{N}\sum_{p=1}^N g_p(t),\; g=x,y,z
```

# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `N`: size of Monte-Carlo sampling, default at 100
- `axis::Int`: , 1,2,3, representing x,y,z axis, set to 3 by default 
- `returnerr::Bool`: wether to return the variance in sampling
"""
function rabi(t::AbstractVector{<:Real}, D::Vector{<:Real}, h::Real; 
    N=100::Int, axis=3::Int, returnerr=false)
    @assert axis in (1,2,3)
    n=length(D)
    f_sum=zeros(length(t)) # sum
    f_var=copy(f_sum) # square sum
    f_sampling=(_rabix,_rabiy,_rabiz)[axis]

    for i in 1:N
        β_p=sum(rand([1,-1],n).*D) 
        f_p=f_sampling(t,β_p,h)
        f_sum+=f_p
        f_var+= i>1 ? (i*f_p-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return (f_sum/N,f_var/(N-1))
    else
        return f_sum/N
    end
end

function rabi(t::AbstractVector{<:Real}, spins::SpinCluster, h::Real; 
    N=100::Int, axis=3::Int, returnerr=false)
    return rabi(t, spins.couplings, h; N=N, axis=axis, returnerr=returnerr)
end

function _rabiz(t::AbstractVector{<:Real},β::Real,h::Real)
    ω=sqrt(h^2+β^2)/2
    cos2=cos.(ω*t).^2
    return (cos2-(cos2.-1)*(β^2-h^2)/(h^2+β^2))/2
end

function _rabiy(t::AbstractVector{<:Real},β::Real,h::Real)
    Ω=sqrt(h^2+β^2)
    return -sin.(Ω*t)/2*h/Ω
end

function _rabix(t::AbstractVector{<:Real},β::Real,h::Real)
    ω=sqrt(h^2+β^2)/2
    sin2=sin.(ω*t).^2
    return sin2*(β*h)/(h^2+β^2)
end


function averagerabi(ensemble::SpinEnsemble, h::Real; M=1000::Int, n_t=200::Int, scale=1.2::Real)
    t=dephasingtime(ensemble, n_t; scale=scale)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=rabi(t, randcoefs(ensemble), h)
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if returnerr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
end
