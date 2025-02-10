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
module MonoAtomic

using UnPack, TOML

export system, initalize!
#This is a structure that holds the data for a molecualr dynamics and Monte Carlo simulation of a monoatomic fluid.
# You initialize it one of three ways:
# 1. system(N=1000, ρ=0.8, T=1.5, dt=0.001) #Lennard-Jones fluid for MD
# 2. system(N=1000, ϕ=0.8) #Hard sphere fluid for MC.
# 3. system("parameters.toml") #Read the parameters from a TOML file.

mutable struct config
    r::Array{Float64, 2}
    v::Array{Float64, 2}
    f::Array{Float64, 2}
end
struct system
    N::Int64
    ρ::Float64
    T::Float64
    L::Float64
    dt::Float64
    atoms::config
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
        atoms_ = config(r_, v_, f_)
        new(N_, ρ_, T_, L_, dt_, atoms_)
    end
    function system(filename::String)
        parameters = TOML.parsefile(filename)
        @unpack N, density, T, dt = parameters
        system(N=N, ρ=density, T=T, dt=dt)
    end
end

function cubic_xtal!(sys::system)
    #Initialize the positions of the particles in a cubic lattice.
end

function init_velocities!(sys::system)
    #Initialize the velocities of the particles from a Gaussian,
    @views begin #This is how you make in-place modifications in Julia.
        v = sys.atoms.v
        v .= randn(Float64, 3, sys.N)
        #Calculate the center of mass velocity
        vcm = sum(v, dims=2) / sys.N
        #Subtract the center of mass velocity from the velocities
        v .-= vcm
        #Scale the velocities to the desired temperature
        v *= sqrt(sys.T)
    end
end

function initalize!(sys::system)
    cubic_xtal!(sys)
    init_velocities!(sys)
end

end