using Random

"""
    rand_locs(N, dim, a)

General method to generate a set of n dimensional locations array. 
a matrix of size (N,d) N random location vectors distributed in a d-dimensional cube, scaled by a at range (-a,a)

# Arguments
- `N`: numer of locations 
- `dim`: dimension 
- `a`: scaling factor
"""
rand_locs(N::Int,dim::Int, a=1.0::Real)=2*a*(rand!(zeros(N,dim)).-1/2)

"""
    randlocscubic(a, b; N=1, dim=3)

Generate `N`` random location vectors distributed in a 3D cube, scaled by a at range `(-b,-a)∪(a,b)``

# Arguments
- `a`: lower bound of sampling range
- `b`: upper bound of sampling range
- `N`: numer of locations, in `[1,2,3]`
- `dim`: dimension of the space
"""
function randlocscubic(a::Real, b::Real; N=1::Int, dim=3)
    @assert dim<=3

    M=zeros(N,3)
    rand!(@view(M[:,1:dim]),[1,-1])
    @view(M[:,1]).*=rand(N).*(b-a).+a
    if dim>1
        @view(M[:,2:dim]).*=rand(N,dim-1).*b
        for i in 1:N
            shuffle!(@view M[i,1:dim])
        end
    end
    M
end


"""
    randlocsspherical(r_min, r_max; N)

Generate `N` random location vectors distributed in a 3D sphere

# Arguments
- `r_min`: lower bound of sampling radius
- `r_max`: upper bound of sampling radius
- `N`: numer of locations 
"""
function randlocsspherical(r_min=0.0::Real, r_max=1.0::Real; N=1::Int)
    M=zeros(N,3)
    r=cbrt.(rand(N).*(r_max^3-r_min^3).+r_min^3)
    ϕ=rand(N).*2pi
    θ=acos.(2*rand(N).-1)

    M[:,1]=r.*sin.(θ).*cos.(ϕ)
    M[:,2]=r.*sin.(θ).*sin.(ϕ)
    M[:,3]=r.*cos.(θ)

    M
end

"""
randlocsspherical(r_min, r_max; N)

Generate `N` random location vectors distributed in a 2D plate

# Arguments
- `r_min`: lower bound of sampling radius
- `r_max`: upper bound of sampling radius
- `N`: numer of locations 
Return:
    a matrix of size (N,3), N random location vectors distributed in a 3D sphere
"""
function randlocspolar(r_min=0.0::Real, r_max=1.0::Real; N=1::Int)
    M=zeros(N,3)
    r=sqrt.(rand(N).*(r_max^2-r_min^2).+r_min^2)
    ϕ=rand(N).*2pi
    
    M[:,1]=r.*cos.(ϕ)
    M[:,2]=r.*sin.(ϕ)      

    M
end