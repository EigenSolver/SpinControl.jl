using DiffEqNoiseProcess

function OrnsteinUhlenbeckNoise(n::Int, dt::Real, params::Union{Tuple,AbstractVector{<:Real}}, X₀::Real = 0)::NoiseProcess
    θ, μ, σ, = params
    if X₀ == -1
        X₀ = σ*randn()
    end
    OU=OrnsteinUhlenbeckProcess(θ,μ,σ,0.0, convert(Float64, X₀))
    OU.dt=dt
    u = nothing; p= nothing;
    calculate_step!(OU, dt, u, p)
    for i in 1:n
        accept_step!(OU,dt, u, p)
    end
    return OU
end


function unitary(noise::NoiseProcess, axis::Int=3)::Vector{<:Matrix}
    n=[0,0,0]
    n[axis]=1
    dt=noise.dt 
    phi_acc=cumsum(noise.u)*dt
    return [rotation(phi, n) for phi in phi_acc]
end

function fid(noise::NoiseProcess)
    dt=noise.dt 
    phi=cumsum(noise.u)*dt
    return cos.(phi)
end

function rabi(noise::NoiseProcess, h::Union{<:Real,AbstractVector{<:Real}})
    """
    Rabi oscillation with noise trace at z direction. Driving field at arbitrary vector h. 
    """
    ρ=[1 0; 0 0];
    f(ρ)=tr(ρ*σ_z)|> real
    
    if h isa Real
        h=[h,0,0]
    end
    
    n=length(noise.u)
    proj_z=zeros(n)
    for i in 1:n 
        proj_z[i]=f(ρ)
        U=rotation(h.+ [0, 0, noise[i]], noise.dt)
        ρ=U*ρ*U'
    end

    return proj_z
end