# SpinEnsembles.jl Documentation


## Functions


```@docs
dipolarcoef(r::AbstractArray{<:Real},z0::AbstractArray{<:Real})
```

```@docs
dipolarcoefs(locs::Matrix{<:Real},z0=[0,0,1.0]::Vector{<:Real})
```

```@docs
randomcoefs(N::Int,dim::Int,a=1::Real)
```

```@docs
fid(t::Real,D::Vector{<:Real})
```

```@docs
averagefid(t::AbstractVector{<:Real}, n_D::Int, sampling_D)
```


```@docs
betasampling(D_set::Vector{<:Real})
```

