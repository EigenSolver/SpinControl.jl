
@doc raw"""
dipolarcoef(r, z0)

Calculate the dipolar interaction strength given by vector `r` and default field at `z0`

```math
D_{ij}=\frac{1-3\cos^2(\theta_{ij})}{2r_{ij}^3} \gamma^2 \hbar
```

# Arguments
- `r`: point vectors from one loc to another    
- `z0`: direction of the background magnetic field 
"""
function dipolarcoef(r::AbstractVector{<:Real},z0::AbstractVector{<:Real})
# suppose z0 is already normalized
cosθ=dot(r,z0)/norm(r) #calculate the cos(θ) between the vector and the z axis
D=0.5*(1-3cosθ^2)/norm(r)^3
return D
end

"""
    dipolarcoefs(locs, z0)

Get a list of dipolar coupling strength between the centered spin and bath

# Arguments
- `locs`: an array of vector, distance from the central spin to the spins in bath 
- `z0`: the direction of external field, set to z axis by default
"""
function dipolarcoefs(locs::Matrix{<:Real},z0=[0,0,1.0]::AbstractVector{<:Real})
    normalize!(z0)
    return map(x->dipolarcoef(x,z0), eachrow(locs))
end

"""
    randcoefs(n, dim, R)
    randcoefs(n, dim, bound; method)

Randomly distribute a spin bath, generate the dipolar coupling strength between the centered spin and bath, 
totally `n` spins are uniformly distributed in a `dim` dimensional cubic space with lenght `a`

# Arguments
- `n`: number of spins in ensemble
- `dim`: dimension
- `R`: scale of ensemble
- `bound::Tuple{Real, Real}`: tuple, indicate bounds of sampling range 

# Options
- `method`: constant, `:spherical` for spherical coordinates or `:cubic` for Cartesian coordinates
"""
function randcoefs(n::Int, dim::Int, bound::Tuple{Real,Real}; method=:cubic)
    return dipolarcoefs(randlocs(n, dim, bound; method=method))
end

randcoefs(n::Int,dim::Int,R=1::Real; method=:cubic)=randcoefs(n,dim,(0,R);method=method)


