#Run this script to generate the data for Problem 3.5.

include("HardSphereMC.jl")
using .HardSphereMC

using Plots

dr = 0.05 #Choose this so that the curve is smooth but small enough to resolve changes in g(r).
sweeps = 10000
equil_time = 2000
calculate_freq = 10

# Create system and RDF calculator
sys = hard_spheres("parameters.toml")
bins = round(Int, box_length(sys)/dr)
calc = rdf(bins=bins)

println("Calculating g(r) for $(number_of_particles(sys)) particles.")
# Run MC simulation
mc_simulation!(
    sys=sys,
    sweeps=sweeps,
    calculations=[calc],
    equil_time=equil_time,
    calculate_freq=calculate_freq
)

# Create plot of g(r)
r = [(i-0.5)*calc.dr for i in 1:calc.bins]
p = plot(r, calc.hist,
    xlabel="r",
    ylabel="g(r)",
    label="Radial Distribution Function",
    linewidth=2)

display(p)
savefig(p, "rdf.png")

# Calculate n(r) using numerical integration
function n(r_values, g_values, dr, ρ)
    n_values = zeros(length(r_values))
    for i in 2:length(r_values)
        # Trapezoidal integration
        y_squared = r_values[1:i].^2
        integrand = y_squared .* g_values[1:i]
        n_values[i] = 4π * ρ * sum(integrand[1:end-1] .+ integrand[2:end]) * dr / 2
    end
    return n_values
end

# Calculate n(r)
ρ = number_of_particles(sys)/box_length(sys)^3
n_values = n(r, calc.hist, calc.dr, ρ)

# Create plot with two y-axes
p2 = plot(r, calc.hist,
    xlabel="r",
    ylabel="g(r)",
    label="g(r)",
    linewidth=2,
    legend=:right)
    
plot!(twinx(), r, n_values,
    ylabel="n(r)",
    label="n(r)",
    color=:red,
    linewidth=2,
    legend=:right,
    legend_position=(0.85, 0.5))

display(p2)
savefig(p2, "rdf_and_n.png")

println("n(L)  ", n_values[end])
# Find first maximum
r1 = findmax(calc.hist)[2]
# Find first minimum after r1
function find_first_local_minimum(data, start_idx)
    for i in start_idx:length(data)-1
        if data[i] < data[i-1] && data[i] < data[i+1]
            return i
        end
    end
    return start_idx
end
r2 = find_first_local_minimum(calc.hist, r1+1)
println("First solvation shell occurs at r = ", r[r2])
r_solv = r[r2]
idx = findmin(abs.(r .- r_solv))[2]
println("Number of neighbors in the first solvation shell ≈ ", n_values[idx])

#write the previous two print statements to file gr_analysis.dat
open("gr_analysis.dat", "w") do file
    write(file, "First solvation shell occurs at r = $(round(r[r2], digits=4))\n")
    write(file, "Number of neighbors in the first solvation shell ≈ $(round(n_values[idx], digits=4))\n")
    write(file, "n(L) = $(n_values[end])\n")
end
