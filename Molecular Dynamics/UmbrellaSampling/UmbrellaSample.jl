includet("UmbrellaSampling.jl")
includet("Potentials.jl")
using .Potentials
using .UmbrellaSampling
# Example usage:
x0, xN = -3.0, 3.0
nwindows = 20
k = 10.0
nsteps = 1e5
Δx = 0.1

umbrellas = umbrella_sample(x0, xN, nwindows, k, quartic, nsteps, Δx, nbins=100)
centers, pmf = wham(umbrellas)

using Plots
plot(centers, pmf, label="PMF", xlabel="x", ylabel="βF(x)")
plot!(centers, quartic.(centers), label="Quartic Potential")