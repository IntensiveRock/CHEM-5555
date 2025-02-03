module vanderPolOscillator

include("ODE_Integrators.jl")
using .ODE_Integrators
using ..ODE_Integrators: ODEIntegrator, solve, ODEProblem
#
#import ..ODE_Integrators: ODEIntegrator, ODEProblem
#import .ODE_Integrators: Heun, RK4
#import .JDE_ODE_Base: solve

using Test, Plots

#Define the van der Pol oscillator problem using the structures (struct) and functions defined in JDE_ODE.jl.

function f(; μ::Float64=1.0)
    return (y,t) -> [y[2], μ*(1 - y[1]^2)*y[2] - y[1]]
end

μ = 1.0
x₀ = 1.0
v₀ = 0.0
Γ = [x₀, v₀]

problem = ODEProblem(f(μ=μ), Γ, 100.0,0.01)

#Solve the problem using RK4
solution = solve(problem, RK2())
plot(solution.t, [y[1] for y in solution.y], label="x(t)", xlabel="t", ylabel="x(t)", title="Van der Pol Oscillator", lw=2)
#write that plot to a file
savefig("vanderPolOscillator.png")


end
