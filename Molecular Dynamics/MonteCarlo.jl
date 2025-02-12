module MonteCarlo

# Perform a Monte Carlo simulation of a fluid. Compute functions from the configuraitons 
# and keep averages.
include("Accessors.jl")
using .Accessors


dot = (x,y) -> sum(x .* y)

function pbc(i::Int64, j::Int64, sys::hard_spheres)
    # Calculate the distance between particles i and j with periodic boundary conditions.
    r = positions(sys)
    L = box_length(sys)
    Linv = 1.0/L
    rᵢⱼ = r[:, i] - r[:, j]
    # Apply periodic boundary conditions
    rᵢⱼ = rᵢⱼ .- L * round.(Int64, rᵢⱼ .* Linv)
    # Compute the squared distance
    dᵢⱼ² = dot(rᵢⱼ, rᵢⱼ)
    return rᵢⱼ, dᵢⱼ²
end

function move!(i::Int64, sys::hard_spheres)
    # Make a trial move on particle i. Accept it or reject it according to Metropolis MC.
    # Compute the distances between particle i and all of its neighbors. If any *Squared* distance is less than dr^2, reject the move.
    # Otherwise, accept it.
    rᵢ = sys.system.atoms.r[:, i]
    accept = 0
    #... Loop over particles j ≠ i, compute rᵢⱼ and dᵢⱼ². If dᵢⱼ² < 1.0, reject the move, otherwise accept it.
    # Set accept = 1 if the move is accepted and 0 otherwise.
    return accept
end

function sweep!(sys::hard_spheres)
    # Perform a Monte Carlo sweep, choosing N particles to move at random and attempting moves for each.
    acceptance_rate = 0.0
    for _ in 1:sys.system.N
        i = rand(1:sys.system.N)
        acceptance_rate += move!(i, sys)
    end
    acceptance_rate /= sys.system.N
    return acceptance_rate
end

function adjust_dr!(sys::hard_spheres, δ::Float64, mean_rate::Float64)
    # Adjust the move size dr to achieve target acceptance rate:
    # If acceptance rate is too low, decrease dr; if too high, increase dr
    if mean_rate <= sys.acceptance
        sys.dr *= (1 + δ)
    else
        sys.dr *= (1 - δ)
    end
end

function mc_simulation!(sys::hard_spheres, sweeps::Int64=10000,f::Function, update_freq::Int64=100, adjust_interval::Int64=50, δ::Float64=0.01)
    # Do a MC simulation of the system, calling the function f every update_freq sweeps and adjusting dr by scaling it by δ.
    # either increasing it or decreasing it multiplicatively.
    accept = 0.0
    for t in 1:sweeps
        accept += sweep!(sys)
        if t % update_freq == 0
            f(sys)
        end
        if t % adjust_interval == 0
            mean_rate = accept/adjust_interval
            @show mean_rate
            adjust_dr!(sys, δ, mean_rate)
            accept = 0.0
        end
    end
end

end