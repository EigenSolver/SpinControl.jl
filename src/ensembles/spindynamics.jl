# common method
function relevanttime(spins::Union{SpinCluster,SpinEnsemble}, n_t::Int; scale::Real = 1.0)
    T2 = coherencetime(spins) * scale
    return 0:T2/n_t:T2
end

# dynamics
"""
fid based on sampling of beta
"""
function fid(
    t::AbstractVector{Float64},
    β::AbstractVector{Float64},   
    h::Real = 0;
    geterr::Bool = false,
)
    N=length(β)
    f_sum = zeros(length(t)) # sum
    f_var = copy(f_sum) # square sum
    @showprogress for i in 1:N
        β_p=β[i]
        ω_p = sqrt(h^2 + β_p^2) / 2
        cos2_p = cos.(ω_p * t) .^ 2
        f_p = (cos2_p + (cos2_p .- 1) * (β_p^2 - h^2) / (h^2 + β_p^2)) / 2
        f_sum += f_p
        f_var += i > 1 ? (i * f_p - f_sum) .^ 2 / (i * (i - 1)) : f_var
    end
    if geterr
        return (f_sum / N, f_var / (N - 1))
    else
        return f_sum / N
    end
end

"""
rabi based on sampling of beta
"""
function rabi(
    t::AbstractVector{Float64},
    β::AbstractVector{Float64}, 
    h::Real;
    axis::Int = 3,
    geterr::Bool = false,
)
    @assert axis in (1, 2, 3)

    N=length(β)
    f_sum = zeros(length(t)) # sum
    f_var = copy(f_sum) # square sum
    f_sampling = (_rabix, _rabiy, _rabiz)[axis]

    @showprogress for i in 1:N
        β_p = β[i]
        f_p = f_sampling(t, β_p, h)
        f_sum += f_p
        f_var += i > 1 ? (i * f_p - f_sum) .^ 2 / (i * (i - 1)) : f_var
    end
    if geterr
        return (f_sum / N, f_var / (N - 1))
    else
        return f_sum / N
    end
end

function rabiperiod(β::AbstractVector{Float64}, h::Real = 0; 
    λ::Real = 0.1, L::Int = 20) # short length fitting 
    
    t0 = π/h
    t = LinRange(t0*(1-λ),t0*(1+λ), L)
    curve=-rabi(t, β, h; axis=2)
    # linear regression
    for i in 2:L
        if curve[i-1]>0 && curve[i]<0
            return t[i-1]+curve[i-1]/(curve[i-1]-curve[i])*(t[i]-t[i-1])
        end
    end
    error("zero points not detected")
end

@doc raw"""
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
function fid(
    t::AbstractVector{Float64},
    ensemble::SpinEnsemble,
    h::Real = 0;
    M::Int = 200,
    N::Int = 100,
    geterr::Bool = false,
)
    β=betasampling(ensemble, M, N)
    return fid(t,β ,h; geterr=geterr)
end

function rabi(
    t::AbstractVector{Float64},
    ensemble::SpinEnsemble,
    h::Real;
    M::Int = 200,
    N::Int = 100,
    axis::Int = 3,
    geterr::Bool = false,
)
    β=betasampling(ensemble, M, N)
    return rabi(t,β,h, axis=axis, geterr=geterr)
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
function fid(
    t::AbstractVector{Float64},
    cluster::SpinCluster,
    h::Real = 0;
    N::Int = 100,
    geterr::Bool = false,
)
    D = cluster.couplings

    if h == 0
        return [mapreduce(cos, *, D * τ) / 2 for τ in t]
    else
        β=betasampling(cluster, N)
        return fid(t,β,h, geterr=geterr)
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
- `geterr::Bool`: wether to return the variance in sampling
"""
function rabi(
    t::AbstractVector{Float64},
    cluster::SpinCluster,
    h::Real;
    N::Int = 100,
    axis::Int = 3,
    geterr::Bool = false,
)
    @assert axis in (1, 2, 3)

    β=betasampling(cluster, N)
    return rabi(t, β, h, axis=axis, geterr=geterr)
end


function _rabiz(t::AbstractVector{<:Real}, β::Real, h::Real)
    Ω = sqrt(h^2 + β^2)
    return (β^2 .+ h^2 * cos.(Ω * t) ) / (2*Ω^2)
end

function _rabiy(t::AbstractVector{<:Real}, β::Real, h::Real)
    Ω = sqrt(h^2 + β^2)
    return -sin.(Ω * t)* h / (2Ω)
end

function _rabix(t::AbstractVector{<:Real}, β::Real, h::Real)
    Ω = sqrt(h^2 + β^2)
    return β*h * (1 .- cos.(Ω * t))  / (2*Ω^2)
end

"""
Get the average driving axis and average driving phase (Rabi frequency) for a spin cluster
""" 
function rabisampling(h::Vector{<:Real}, cluster::SpinCluster; N::Int=100)
    z0=cluster.ensemble.z0
    β = betasampling(cluster, N)
    h_p= β.*z0' .+ h'
    
    return h_p
end

function rabisampling(h::Real, cluster::SpinCluster, 
    aim::Vector{<:Real}=[1,0,0]; N::Int=100)
    normalize!(aim)
    h_p = rabisampling(h*aim, cluster, N=N)
    
    Ω_p = sqrt.(sum(abs2, h_p, dims=2))
    return vec(Ω_p), h_p./Ω_p
end

function rabisampling(h::Real, t::Real, cluster::SpinCluster, 
    aim::Vector{<:Real}=[1,0,0]; N::Int=100) 
    Ω_p, n_p = rabisampling(h, cluster, aim, N=N)
    return t*Ω_p, n_p
end 

"""
Find the Rabi period of the given ensemble under driving field, using linear regression at slope.
"""
function rabiperiod(ensemble::SpinEnsemble, h::Real = 0; 
    M::Int = 1000, N::Int = 100, λ::Real = 0.1, L::Int = 20) # short length fitting 
    β=betasampling(ensemble, M, N)
    return rabiperiod(β,h,λ=λ,L=L)
end