

"""
Args:
    N: number of point generated
    dim: dimension of the space
    a: scaling factor 
Return:
    an array of N random location vectors, in given dimension, scaled by a at range [-a,a]
"""
function rand_locs_cubic(a::Real, b::Real; projection=:false)
    tmp=rand([-1,1],3).*rand(Float64,3).*(b-a).+a
    if projection
        tmp[end]=0.0
    end
    tmp
end

function rand_locs_spherical(r_max=1.0::Real, r_min=0.0::Real; N=1::Int, projection=:false)
    r=rand(N).*(r_max-r_min).+r_min
    ϕ=rand(N).*2pi
    if projection
        θ=ones(N).*pi/2
    else
        θ=rand(N).*pi
    end
    collect(zip(r.*sin.(θ).*cos.(ϕ),r.*sin.(θ).*sin.(ϕ),r.*cos.(θ)))
end

function _rand_loc_spherical(r_max=1.0::Real, r_min=0.0::Real; projection=:false)
    r=rand()*(r_max-r_min)+r_min
    ϕ=rand()*2pi
    if projection
        θ=pi/2
    else
        θ=rand()*pi
    end
    [r*sin(θ)*cos(ϕ),r*sin(θ)*sin(ϕ),r*cos(θ)]
end