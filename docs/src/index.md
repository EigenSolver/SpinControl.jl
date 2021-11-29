# SpinEnsembles.jl Documentation


## Functions

```@docs
bath_vectors(loc0::Vector{<:Real},loc_bath::Matrix{<:Real})=map(x->x-loc0,eachrow(loc_bath))
```

```@docs
function dipolar_coef(r::AbstractArray{<:Real},z0::AbstractArray{<:Real})
```

<!-- ## Index

```@index
``` -->