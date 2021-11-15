using Random

"""
General method to generate a set of n dimensional locations array
Args:
    N: numer of locations 
    dim: dimension 
    a: scaling factor
Return:
   a matrix of size (N,d) N random location vectors distributed in a d-dimensional cube, scaled by a at range (-a,a)
"""
rand_locs(N::Int,dim::Int, a=1.0::Real)=2*a*(rand!(zeros(N,dim)).-1/2)

"""
Args:
    a: lower bound of sampling range
    b: upper bound of sampling range
    N: numer of locations 
    projection: a boolean to decide whether to project the points to a 2D plane
Return:
   a matrix of size (N,3) N random location vectors distributed in a 3D cube, scaled by a at range (-b,-a)∪(a,b)
"""
function rand_locs_cubic(a::Real, b::Real; N=1::Int, dim=3)
    M=zeros(N,3)
    rand!(@view(M[:,1:dim]),[1,-1])
    M.*=rand(N,3).*(b-a).+a
    M
end


"""
Args:
    r_min: lower bound of sampling radius
    r_max: upper bound of sampling radius
    N: numer of locations 
    projection: a boolean to decide whether to project the points to a 2D plane
Return:
    a matrix of size (N,3), N random location vectors distributed in a 3D sphere
"""
function rand_locs_spherical(r_max=1.0::Real, r_min=0.0::Real; N=1::Int, projection=:false)
    M=zeros(N,3)
    r=rand(N).*(r_max-r_min).+r_min
    ϕ=rand(N).*2pi
    if projection
        θ=ones(N).*pi/2
    else
        θ=rand(N).*pi
    end
    
    M[:,1]=r.*sin.(θ).*cos.(ϕ)
    M[:,2]=r.*sin.(θ).*sin.(ϕ)
    M[:,3]=r.*cos.(θ)

    M
end