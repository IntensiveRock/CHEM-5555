module MonoAtomic
"""
Monoatomic module that provides core functionality for molecular dynamics and Monte Carlo 
simulations of a monoatomic fluid.

system

Represents a molecular dynamics system with positions, velocities, and forces.

Fields:
- N::Int64: Number of particles
- ρ::Float64: Number Density
- T::Float64: Temperature
- L::Float64: Box length
- r::Array{Float64, 2}: Positions (3×N)
- v::Array{Float64, 2}: Velocities (3×N)
- f::Array{Float64, 2}: Forces (3×N)
- dt::Float64: Time step


"""
using UnPack, TOML

export system, initalize!

function density(packing_fraction::Float64)
    return 6*packing_fraction/π
end

function packing_fraction(density::Float64)
    return π*density/6
end

#This is a structure that holds the data for a molecualr dynamics and Monte Carlo simulation of a monoatomic fluid.
# You initialize it one of three ways:
# 1. system(N=1000, ρ=0.8, T=1.5, dt=0.001) #Lennard-Jones fluid for MD
# 2. system(N=1000, ϕ=0.8) #Hard sphere fluid for MC.
# 3. system("parameters.toml") #Read the parameters from a TOML file.

mutable struct system
    N::Int64
    ρ::Float64
    T::Float64
    L::Float64
    r::Array{Float64, 2}
    v::Array{Float64, 2}
    f::Array{Float64, 2}
    dt::Float64
    function system(; N::Int64, ϕ::Float64)
        N_ = N
        ρ_ = density(ϕ)
        L_ = (N_/ρ_)^(1/3)
        T_ = 0.0
        r_ = zeros(Float64, 3, N_)
        v_ = zeros(Float64, 3, N_)
        f_ = zeros(Float64, 3, N_)
        new(N_, ρ_, T_, L_, r_, v_, f_, 0.0)
    end
    function system(; N::Int64=1000, ρ::Float64=0.8, T::Float64=1.5, dt::Float64=0.001)
        N_ = N
        ρ_ = ρ
        T_ = T
        dt_ = dt
        L_ = (N_/ρ_)^(1/3)
        r_ = zeros(Float64, 3, N_)
        v_ = zeros(Float64, 3, N_)
        f_ = zeros(Float64, 3, N_)
        dt_ = 0.001
        new(N_, ρ_, T_, L_, r_, v_, f_, dt_)
    end
    function system(filename::String)
        parameters = TOML.parsefile(filename)
        @unpack N, density, T, dt = parameters
        if T == 0.0
            return system(N=N, ϕ = density)
        end
        ρ_ = density
        L_ = (N / ρ_)^(1 / 3)
        r_ = zeros(Float64, 3, N_)
        v_ = zeros(Float64, 3, N_)
        f_ = zeros(Float64, 3, N_)
        new(N, ρ_, T, L_, r_, v_, f_, dt)
    end
end

function cubic_xtal!(sys::system)
    #Initialize the positions of the particles in a cubic lattice.
end

function init_velocities!(sys::system)
    #Initialize the velocities of the particles from a Gaussian,
    for i in 1:sys.N
        sys.v[1, i] = randn()
        sys.v[2, i] = randn()
        sys.v[3, i] = randn()
    end
    #Calculate the center of mass velocity
    vcm = zeros(Float64, 3)
    for i in 1:sys.N
        vcm += sys.v[:, i]
    end
    vcm /= sys.N
    #Subtract the center of mass velocity from the velocities
    for i in 1:sys.N
        sys.v[:, i] -= vcm
    end
    #Scale the velocities to the desired temperature
    sys.v*=sqrt(sys.T)
end

function initalize!(sys::system)
    cubic_xtal!(sys)
    init_velocities!(sys)
end

end