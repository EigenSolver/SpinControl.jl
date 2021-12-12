# Functions

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
randcoefs(n::Int, dim::Int, bound::Tuple{Real,Real}; method=:cubic)
```

```@docs
fid(t::Real,D::Vector{<:Real})
```

```@docs
fid(t::AbstractVector{<:Real},D::Vector{<:Real},h::Real; N=100::Int, geterr=:false)
```

```@docs
averagefid(t::AbstractVector{<:Real}, M::Int, sampling_D::Function)
```

```@docs
rabi(t::AbstractVector{<:Real}, D::Vector{<:Real}, h::Real; N=100::Int, axis=3::Int, returnerr=:false)
```

```@docs
rabi(t::AbstractArray, h::Real, b::Real; axis=3)
```

```@docs
betasampling(D_set::Vector{<:Real})
```

```@docs
decaytime(D::Vector{<:Real},M::Int=500;len=500::Int,n_sigma=2::Real)
```
