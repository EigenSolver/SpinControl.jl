# dynamics
@doc raw"""
    fid(t, M, sampling_D; [options]...)
    fid(t, M, sampling_D; [options]...)

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
function fid(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h=0::Real; 
    M=200::Int, N=100::Int, geterr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    cluster=SpinCluster(ensemble)
    @showprogress for i in 1:M
        f_d=fid(t, cluster, h; N=N)
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
        reroll!(cluster)
    end
    if geterr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function fid(ensemble::SpinEnsemble, h=0::Real;
    M=200::Int, N=100::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=relevanttime(ensemble, n_t; scale=scale)
    return fid(t, ensemble, h; M=M, N=N, geterr=geterr)
end

function rabi(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h::Real;
    M=200::Int, N=100::Int, axis=3::Int, geterr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    cluster=SpinCluster(ensemble)
    @showprogress for i in 1:M
        f_d=rabi(t, cluster, h; N=N, axis=axis)
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
        reroll!(cluster)
    end
    if geterr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function rabi(ensemble::SpinEnsemble, h::Real; 
    M=200::Int, n_t=200::Int, scale=1.0::Real, N=100::Int, axis=3::Int, geterr=false)
    T2=coherencetime(spins)*scale
    return rabi(t, ensemble, h; M=M, N=N, axis=axis, geterr=geterr)
end

# common method
function relevanttime(spins::Union{SpinCluster,SpinEnsemble}, n_t::Int; scale=1.0::Real)
    T2=coherencetime(spins)*scale
    return 0:T2/n_t:T2 
end

@doc raw"""
    fid(t, D, h; N)

Calculate the free induction dacay of the given spin cluster,
with given transverse magnetic field, using Monte-Carlo sampling

```math
f_p(t)=\frac{1}{2}[\cos^2(\omega_p t)+\sin^2(\omega_p t) (n_x^2-n_z^2)]
```

```math
\bar{f}(t)=\frac{1}{N}\sum_{p=1}^N f_p(t)
```

when h=0, the equation is reduced to 

```math
f(t)=\frac{1}{2}\prod_j \cos(D_jt)
```

# Arguments
- `t`: discrete array marking time 
- `D`: a set of the coupling strengths
- `h`: strength of transverse field 

# Options
- `N`: number of Monte-Carlo sampling
- `geterr::Bool`: wetehr to return the error of the monte-sampling 
"""
function fid(t::AbstractVector{Float64}, cluster::SpinCluster,h=0::Real; 
    N=100::Int, geterr=false)
    D=cluster.couplings

    if h==0
        return [mapreduce(cos,*,D*τ)/2 for τ in t]
    end

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

function fid(cluster::SpinCluster,h=0::Real; 
    N=100::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=relevanttime(cluster,n_t; scale=scale)
    fid(t, cluster, h; N=N, geterr=geterr)
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
- `geterr::Bool`: wether to return the variance in sampling
"""
function rabi(t::AbstractVector{Float64}, cluster::SpinCluster, h::Real; 
    N=100::Int, axis=3::Int, geterr=false)
    @assert axis in (1,2,3)
    
    D=cluster.couplings
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
    if geterr
        return (f_sum/N,f_var/(N-1))
    else
        return f_sum/N
    end
end

function rabi(cluster::SpinCluster, h::Real; 
    N=100::Int, axis=3::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=relevanttime(cluster,n_t; scale=scale)
    @assert (t[2]-t[1])<π/(20*h)
    return rabi(t, cluster, h; N=N, axis=axis, geterr=geterr)    
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
