# SpinEnsembles.jl

## Functions

```@docs
dipolarcoef(r::AbstractVector{<:Real},z0::AbstractVector{<:Real})
```

```@docs
dipolarcoefs(locs::Matrix{<:Real},z0=[0,0,1.0]::Vector{<:Real})
```

```@docs
dipolarlinewidth(D::Vector{<:Real})
```

```@docs
randcoefs(N::Int,dim::Int,a=1::Real)
```

```@docs
fid(t::Real,D::Vector{<:Real})
```

```@docs
fid(t::Real,D::Vector{<:Real},h::Real;N=1::Int)
```

```@docs
averagefid(t::AbstractVector{<:Real}, n_D::Int, sampling_D)
```

```@docs
betasampling(D_set::Vector{<:Real})
```

```@docs
decaytime(D::Vector{<:Real},M::Int=500;len=500::Int,n_sigma=2::Real)
```
