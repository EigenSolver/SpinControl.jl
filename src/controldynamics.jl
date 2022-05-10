import ProgressMeter: @showprogress

function unitary(pulse::SquarePulse, β::Real=0, z0::Vector{<:Real}=[0,0,1])
    return rotation(pulse.aim.*pulse.h + z0.*β, pulse.t)
end

function unitary(pulse::Idle, β::Real=0, z0::Vector{<:Real}=[0,0,1])
    return rotation(z0.*β, pulse.t)
end


function krausoperators(pulse::Pulse, β::AbstractVector{<:Real},
    c::Vector{<:Real}=normalize!(ones(size(β)), 1),
    z0::Vector{<:Real}=[0,0,1]
    )::Vector{<:Matrix}
    return sqrt.(c).*[unitary(pulse, β_k, z0) for β_k in β]
end


function unitary(seq::Sequence, β::Real=0, z0::Vector{<:Real}=[0,0,1])
    U0=unitary(seq.idle,β,z0)
    Un=[unitary(g,β,z0) for g in seq.gates]
    V=σ_i
    for i in seq.order
        U= i==0 ? U0 : (sign(i)>0 ? Un[abs(i)] : Un[abs(i)]')
        V=U*V
    end
    return V
end

function krausoperators(seq::Sequence, β::AbstractVector{<:Real},
    c::Vector{<:Real}=normalize!(ones(size(β)), 1),
    z0::Vector{<:Real}=[0,0,1]
    )::Vector{<:Matrix}
    
    return sqrt.(c).*[unitary(seq, β_k, z0) for β_k in β]
end

function deploy(ψ::Union{Vector{ComplexF64},Matrix{ComplexF64}}, seq::Sequence, n::Int, β::Real, 
    z0::Vector{<:Real}=[0,0,1]; cycle::Int=1, gettime::Bool=false)
    dt=seq.idle.t/n
    U0=unitary(Idle(dt),β,z0)
    Up=map(g->unitary(g,β,z0), seq.gates)
    Un=map(g->unitary(g,β,-z0), seq.gates)
    
    t_cycle= cycleslice(seq,n)
    N=length(t_cycle)   
    
    T=typeof(ψ)
    ψ_arr = Vector{T}(undef,1+(N-1)*cycle)
    p=1;
    ψ_arr[p]=ψ

    for k in 1:cycle
        for i in seq.order
            if i==0
                for j in 1:n
                    ψ=evolve(ψ, U0)
                    p+=1
                    ψ_arr[p] = ψ
                end
            else
                U= sign(i)>0 ? Up[abs(i)] : Un[abs(i)]'
                ψ=evolve(ψ, U)
                p+=1
                ψ_arr[p] = ψ
            end
        end
    end
    if gettime 
        t_c=cycletime(seq)
        t_arr = append!(t_cycle, [t_cycle[2:end] .+ i*t_c for i in 1:cycle-1]...)
        return t_arr, ψ_arr
    else
        return ψ_arr
    end
end


function deploy(ρ::Matrix{ComplexF64}, seq::Sequence, n::Int, β::Vector{<:Real}, 
    c::Vector{<:Real}=normalize!(ones(size(β)), 1), z0::Vector{<:Real}=[0,0,1]; cycle::Int=1)
    
    t_cycle= cycleslice(seq,n)
    t_c=cycletime(seq)
    N=length(t_cycle)
    
    ρ_arr = [[0.0im 0; 0 0] for i in 1:(1+(N-1)*cycle)]
    @showprogress for k in 1:length(β)
        ρ_arr= ρ_arr .+ c[k].*deploy(ρ, seq, n, β[k], z0, cycle=cycle, gettime=false)
    end
    
    t_arr = append!(t_cycle, [t_cycle[2:end] .+ i*t_c for i in 1:cycle-1]...)
    return t_arr, ρ_arr
end