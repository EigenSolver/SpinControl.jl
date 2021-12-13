using Random

"""
    randcartesianlocs(n, a, b; dim=3)

Generate `n` random location vectors distributed in a 3D cube, scaled by a at range `(-b,-a)∪(a,b)`

# Arguments
- `a`: lower bound of sampling range
- `b`: upper bound of sampling range
- `n`: numer of locations, in `[1,2,3]`

# Options
- `dim`: dimension of the space
"""
function randcartesianlocs(n::Int, a::Real, b::Real; dim=3::Int)
    @assert dim<=3

    M=zeros(n,3)
    rand!(@view(M[:,1:dim]),[1,-1])
    @view(M[:,1]).*=rand(n).*(b-a).+a
    if dim>1
        @view(M[:,2:dim]).*=rand(n,dim-1).*b
        for i in 1:n
            shuffle!(@view M[i,1:dim])
        end
    end
    M
end

function randcartesianlocs(n::Int, R::Real; dim=3::Int)
    @assert dim<=3
    M=zeros(n,3)
    rand!(@view(M[:,1:dim]))
    @view(M[:,1:dim]).-=1/2
    return 2R*M
end

"""
    randsphericallocs(n, r_min, r_max)

Generate `n` random location vectors distributed in a 3D sphere,
return a matrix of size (n,3), n random location vectors distributed in a 3D sphere

# Arguments
- `r_min`: lower bound of sampling radius
- `r_max`: upper bound of sampling radius
- `n`: numer of locations 
"""
function randsphericallocs(n::Int, r_min=0.0::Real, r_max=1.0::Real)
    M=zeros(n,3)
    r=cbrt.(rand(n).*(r_max^3-r_min^3).+r_min^3)
    ϕ=rand(n).*2pi
    θ=acos.(2*rand(n).-1)

    M[:,1]=r.*sin.(θ).*cos.(ϕ)
    M[:,2]=r.*sin.(θ).*sin.(ϕ)
    M[:,3]=r.*cos.(θ)

    M
end

"""
randsphericallocs(n, r_min, r_max)

Generate `n` random location vectors distributed in a 2D plate,
return a matrix of size (n,3), n random location vectors distributed in a 3D sphere

# Arguments
- `n`: numer of locations 
- `r_min`: lower bound of sampling radius
- `r_max`: upper bound of sampling radius
"""
function randpolarlocs(n::Int, r_min=0.0::Real, r_max=1.0::Real)
    M=zeros(n,3)
    r=sqrt.(rand(n).*(r_max^2-r_min^2).+r_min^2)
    ϕ=rand(n).*2pi
    
    M[:,1]=r.*cos.(ϕ)
    M[:,2]=r.*sin.(ϕ)      

    M
end