module UmbrellaSampling

include("Potentials.jl")
include("Sampling.jl")
using StatsBase: Histogram, fit
using .Potentials, .Sampling

export umbrella, umbrella_sample, wham

mutable struct umbrella
    βΦ::Function       # Energy function
    βw::Function      # Bias function
    edges::Vector{Float64}
    weighted_histogram::Vector{Float64}
    unweighted_histogram::Vector{Float64}
end

# Constructor for umbrella with harmonic bias
function umbrella(βΦ::Function, x0::Float64, k::Float64, edges::Vector{Float64})
    βw(x) = 0.5 * k * (x - x0)^2
    weighted_hist = zeros(length(edges)-1)
    unweighted_hist = zeros(length(edges)-1)
    return umbrella(βΦ, βw, edges, weighted_hist, unweighted_hist)
end

function sample!(u::umbrella, nsteps::Int64, Δx::Float64)
    _, traj = MC(u.edges[1], 
                x -> u.βΦ(x) + u.βw(x),
                nsteps,
                Δx,
                criterion=(ΔβΦ -> metropolis(ΔβΦ)))
    
    h = fit(Histogram, traj, u.edges)
    u.weighted_histogram = h.weights ./ sum(h.weights)
    u.unweighted_histogram = u.weighted_histogram .* exp.(u.βw.([u.edges[i] + 0.5*(u.edges[i+1] - u.edges[i]) for i in 1:length(u.edges)-1]))
end

function umbrella_sample(x0::Float64, xN::Float64, nwindows::Int64, k::Float64, βΦ::Function,
                        nsteps::Int64, Δx::Float64; nbins::Int64=100)
    edges = range(x0, xN, length=nbins+1) |> collect
    windows = range(x0, xN, length=nwindows)
    umbrellas = [umbrella(βΦ, x_center, k, edges) for x_center in windows]
    
    for u in umbrellas
        sample!(u, nsteps, Δx)
    end
    
    return umbrellas
end

function wham(umbrellas::Vector{umbrella}; max_iter::Int64=1000, tol::Float64=1e-6)
    edges = umbrellas[1].edges
    centers = [(edges[i] + edges[i+1])/2 for i in 1:length(edges)-1]
    nwindows = length(umbrellas)
    
    # Initialize free energies
    F = zeros(nwindows)
    P = zeros(length(centers))
    
    # WHAM iteration
    for iter in 1:max_iter
        F_old = copy(F)
        
        # Calculate P(x)
        for j in 1:length(centers)
            # Raw histogram counts from each window
            numerator = sum(u.weighted_histogram[j] for u in umbrellas)
            
            # Proper bias removal with correct sign convention
            denominator = sum(exp(F[i]) * exp(-u.βw(centers[j])) 
                            for (i,u) in enumerate(umbrellas))
            
            P[j] = numerator / denominator
        end
        
        # Update free energies (note sign change)
        for (i,u) in enumerate(umbrellas)
            F[i] = -log(sum(P[j] * exp(-u.βw(centers[j])) 
                          for j in 1:length(centers)))
        end
        
        if maximum(abs.(F - F_old)) < tol
            break
        end
    end
    
    # Final PMF calculation
    pmf = -log.(P)
    pmf .-= minimum(pmf)  # Shift minimum to zero
    
    return centers, pmf
end

end
