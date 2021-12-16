
"""
    SpinEnsemble(n, dim, z0, r, R)
    SpinEnsemble(rho, dim, z0, r, R)

# Arguments:
- `n::Int64`: number of spins in the ensemble 
- `rho::Float64`: density of spins
- `dim::Int64`: dimension of the ensemble 
- `z0::Vector{Float64}`: direction of background field 
- `r::Float64`: minimum radius of distribution area 
- `R::Float64`: maximum radius of distribution area
"""
struct SpinEnsemble
    n::Int
    rho::Float64
    dim::Int
    z0::Vector{Float64}
    r::Float64
    R::Float64
    shape::Symbol

    function SpinEnsemble(n_or_rho::Union{Int,Float64},dim::Int,z0::AbstractVector{<:Real},
        r::Real, R::Real, shape=:spherical::Symbol)
        @assert n_or_rho>0
        @assert dim in (1,2,3)
        @assert shape in (:spherical, :cubic)
        V=_volume(dim, r, R, shape)
        if n_or_rho isa Int
            n=n_or_rho
            rho=n/V
        elseif n_or_rho isa Float64
            rho=n_or_rho
            n=floor(Int, rho*V)
        end
        return new(n,rho,dim,normalize(z0),float(r),float(R),shape)  
    end
end

function _volume(dim::Int, r::Real, R::Real, shape::Symbol)
    if shape==:spherical
        V=((2*(R-r)),(π*(R^2-r^2)),(4π/3*(R^3-r^3)))[dim]
    elseif shape==:cubic
        V=((2*(R-r)),(4*(R^2-r^2)),(8*(R^3-r^3)))[dim]
    end
    return V
end

function volume(ensemble::SpinEnsemble)
    return _volume(ensemble.dim, ensemble.r, ensemble.R, ensemble.shape)
end


# methods & properties
isdilute(ensemble::SpinEnsemble)=ensemble.rho<1
randlocs(ensemble::SpinEnsemble)=randlocs(ensemble.n, ensemble.dim, (ensemble.r,ensemble.R);
 method=ensemble.shape)
randcoefs(ensemble::SpinEnsemble)=dipolarcoefs(randlocs(ensemble))


function coherencetime(ensemble::SpinEnsemble)
    if isdilute(ensemble)
        if ensemble.dim==3
            return 6/(ensemble.rho*π^2)*(3*sqrt(3))/8
        elseif ensemble.dim==2
            return 2*(ensemble.rho*π/2*gamma(1/3))^(-3/2)
        elseif ensemble.dim==1
            return -2*(ensemble.rho/sqrt(3)*gamma(-1/3))^(-3)
        else
            throw(DomainError(ensemble.dim,"Invalide ensemble dimension."))
        end
    else
        throw(DomainError("Ensemble is non-dilute"))
    end
end

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
function fid(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h=0; M=1000::Int, geterr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=fid(t, dipolarcoefs(randlocs(ensemble), h))
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if geterr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function fid(ensemble::SpinEnsemble;
    M=1000::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=relevanttime(ensemble, n_t; scale=scale)
    return fid(t, ensemble; M=M, geterr=geterr)
end

function rabi(t::AbstractVector{Float64}, ensemble::SpinEnsemble, h::Real;
    M=1000::Int, geterr=false)
    f_sum=zeros(length(t))
    f_var=copy(f_sum)
    @showprogress for i in 1:M
        f_d=rabi(t, randcoefs(ensemble), h)
        f_sum+=f_d
        f_var+= i>1 ? (i*f_d-f_sum).^2/(i*(i-1)) : f_var
    end
    if geterr
        return f_sum/M, f_var/(M-1)
    else
        return f_sum/M
    end
end

function rabi(ensemble::SpinEnsemble, h::Real; 
    M=1000::Int, n_t=200::Int, scale=1.0::Real, geterr=false)
    t=relevanttime(ensemble, n_t; scale=scale)
    return rabi(t, ensemble, h; M=M, geterr=geterr)
end

# common method
function relevanttime(spins::Union{SpinCluster,SpinEnsemble}, n_t::Int; scale=1.0::Real)
    T2=coherencetime(spins)*scale
    return 0:T2/n_t:T2 
end

mutable struct SpinCluster
    ensemble::SpinEnsemble
    locations::Matrix{Float64}
    couplings::Vector{Float64}
    linewidth::Float64

    function SpinCluster(ensemble::SpinEnsemble)
        locs=randlocs(ensemble)
        D=dipolarcoefs(locs)
        new(ensemble, locs, D, dipolarlinewidth(D))
    end 

    SpinCluster(locations::Matrix{Float64})=
    (D=dipolarcoefs(locations); new(missing, locations, D, dipolarlinewidth(D)))
end

isdilute(cluster::SpinCluster)=cluster.ensemble.rho<1

@doc raw"""
    betasampling(cluster; N)

Get the Gaussian linewidth of dipolar coupling for the given spin cluster

    ```math
\beta_p = \sum_j p_j \,D_j,\; p_j=\pm 1
```

# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
- `ensemble`: a `SpinEnsemble`, use this function as the method of the type 
# Options
- `N`: size of Monte-Carlo sampling
"""
function betasampling(cluster::SpinCluster; N=1::Int)
    return [sum(rand([1,-1],n).*cluster.couplings) for i in 1:N]
end



@doc raw"""
    dipolarlinewidth(cluster)

Get the linewidth of D, which follows Gaussian distribution. 

```math
b=\sqrt{\sum_j D_j^2}
```
# Arguments
- `D`: coupling strength of dipolar interactions, a array of floats
"""
dipolarlinewidth(cluster::SpinCluster)=sqrt(mapreduce(abs2,+,cluster.couplings)) # in time computed and stored


@doc raw"""
    coherencetime(D, n_t=500; scale=1.0)

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
coherencetime(cluster::SpinCluster)=π/dipolarlinewidth(cluster)

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
