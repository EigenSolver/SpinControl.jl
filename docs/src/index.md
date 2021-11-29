# SpinEnsembles.jl Documentation


## Functions


```@docs
dipolar_coef(r::AbstractArray{<:Real},z0::AbstractArray{<:Real})
```

```@docs
bath_dipolar_coefs(vec_bath::Matrix{<:Real},z0=[0,0,1.0]::Vector{<:Real})
```

```@docs
rand_bath_dipolar_coefs(N::Int,dim::Int,a=1::Real)
```

```@docs
ensemble_FID(t::Real,D::Vector{<:Real})
```

```@docs
beta_sampling(D_set::Vector{<:Real})
```