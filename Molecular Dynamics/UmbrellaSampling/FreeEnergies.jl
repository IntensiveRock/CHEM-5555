using Revise
using Plots
using StatsBase
includet("Potentials.jl")
includet("Sampling.jl")

plotly()

using .Potentials, .Sampling

function free_energy(x, bins_::Int64=100)
    # Create histogram edges
    edges = range(minimum(x), maximum(x), length=bins_+1)
    
    # Fit histogram without normalize keyword
    P = fit(Histogram, x, edges)
    
    # Normalize manually
    weights = P.weights / sum(P.weights)
    
    # Calculate free energy
    βF = -log.(weights)
    
    return edges[1:end-1], βF
end


println("Find the stepsize Δx")
Δx = find_Δx(-1.0, quartic, 1000, criterion=metropolis,acceptance_target=0.5, window_size=100)

γ, x = MC(-1.0, quartic, Int(1e8), Δx, criterion=metropolis)
println("Acceptance rate: ", γ)
println("Optimal Δx: ", Δx)
println("Plotting trajectory starting at -1.0")
histogram(x, bins=50, normalize=true,
    title="Histogram of x values", xlabel="x", ylabel="Frequency")


q, βF = free_energy(x)

plot(q, βF,
    label="Free Energy",
    xlabel="x",
    ylabel="βF(x)")


# Evaluate quartic potential on same grid as free energy
V = quartic.(q)

# Inf values in βF do occur because the histogram can have zero values.
# Find shift by minimizing squared difference between βF and V, excluding Inf values
valid_indices = .!isinf.(βF)
shift = mean((βF .- V)[valid_indices])

# Plot both curves
plot(q, βF .- shift, label="Free Energy", xlabel="x", ylabel="Energy")
plot!(q, V, label="Quartic Potential")