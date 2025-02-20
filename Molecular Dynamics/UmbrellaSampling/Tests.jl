using Test
using Plots
using Revise
includet("Potentials.jl")
using .Potentials
includet("Sampling.jl")
using .Sampling

@testset "Quartic Potential Tests" begin
    # Test find_minimum on quartic potential
    min_x = Potentials.find_minimum(quartic)
    @test isapprox(quartic(min_x), -10.0, atol=1e-6)
    
    # Print minimum
    println("Minimum of quartic potential at x = ", min_x)
    println("Value at minimum = ", quartic(min_x))
end


