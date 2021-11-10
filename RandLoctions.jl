"""
Args:
    a: lower bound of sampling range
    b: upper bound of sampling range
    N: numer of locations 
    projection: a boolean to decide whether to project the points to a 2D plane
Return:
    an array of N random location vectors, in given dimension, scaled by a at range (-b,-a)∪(a,b)
"""
function rand_locs_cubic(a::Real, b::Real; N=1, projection=:false)
    M=zeros(N,3)
    rand!(M,[1,-1])
    M.*=rand(N,3).*(b-a).+a
    if projection
        M[:,end].=0.0
    end
    M
end

