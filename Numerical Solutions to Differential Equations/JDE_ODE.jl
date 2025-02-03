module JDE_ODE

#This is an implementation for several integrators. I use similar structures to those we've been using,
#but now I'm reaching into my bag of tricks and pulling out template functions. This is a way to write
#functions that can be used for different types of arguments. In our case we don't want to have to re-write
#code that we've built and debugged simply because we have more than one variable.

include("JDE_ODE_Base.jl")
using .JDE_ODE_Base
import Roots: fzero
import UnPack: @unpack
import .JDE_ODE_Base: step #This is necessary to import and extend the function step.
using Test

export ForwardEuler, BackwardEuler, Heun, Midpoint, RK2, RK4, ImplicitMidpoint, BackwardEuler, CrankNicolson, Trapezoidal

struct ForwardEuler <: ODEIntegrator end
function step(problem::ODEProblem, ::ForwardEuler, yₙ::Float64, tₙ::Float64)
    @unpack f, h = problem
    return yₙ + h*f(yₙ,tₙ)
end

#Second-order explicit integrators
struct RK2 <: ODEIntegrator end
struct Heun <: ODEIntegrator end
function step(problem::ODEProblem, ::Union{Heun,RK2}, yₙ::typename, tₙ::Float64) where typename
    @unpack f, h = problem
    k₁ = f(yₙ,tₙ)
    k₂ = f(yₙ + h*k₁, tₙ + h)
    return yₙ + (h/2)*(k₁ + k₂)
end

struct Midpoint <: ODEIntegrator end
function step(problem::ODEProblem, ::Midpoint, yₙ::typename, tₙ::Float64) where typename
    @unpack f, h = problem
    k₁ = f(yₙ,tₙ)
    k₂ = f(yₙ + (h/2)*k₁, tₙ + h/2)
    return yₙ + h*k₂
end

#Fourh-order RK4.
struct RK4 <: ODEIntegrator end
function step(problem::ODEProblem, ::RK4, yₙ::typename, tₙ::Float64) where typename
    @unpack f, h = problem
    k₁ = f(yₙ,tₙ)
    k₂ = f(yₙ + (h/2)*k₁, tₙ + h/2)
    k₃ = f(yₙ + (h/2)*k₂, tₙ + h/2)
    k₄ = f(yₙ + h*k₃, tₙ + h)
    return yₙ + (h/6)*(k₁ + 2k₂ + 2k₃ + k₄)
end

#Implicit integrators, only implemented for scalar problems
struct ImplicitMidpoint <: ODEIntegrator end
function step(problem::ODEProblem, ::ImplicitMidpoint, yₙ::Float64, tₙ::Float64)
    @unpack f, h = problem
    g(u) = u - yₙ - h*f((u + yₙ)/2, tₙ + h/2)
    return fzero(g, yₙ + h*f(yₙ,tₙ))
end

struct BackwardEuler <: ODEIntegrator end
function step(problem::ODEProblem, ::BackwardEuler, yₙ::Float64, tₙ::Float64)
    @unpack f, h = problem
    g(u) = u - yₙ - h*f(u,tₙ+h)
    return fzero(g, yₙ + h*f(yₙ,tₙ))
end

struct CrankNicolson <: ODEIntegrator end
function step(problem::ODEProblem, ::CrankNicolson, yₙ::Float64, tₙ::Float64)
    @unpack f, h = problem
    g(u) = u - yₙ - h/2*(f(yₙ,tₙ) + f(u,tₙ+h))
    return fzero(g, yₙ + h*f(yₙ,tₙ))
end

struct Trapezoidal <: ODEIntegrator end
function step(problem::ODEProblem, ::Trapezoidal, yₙ::Float64, tₙ::Float64)
    @unpack f, h = problem
    g(u) = u - yₙ - h/2*(f(yₙ,tₙ) + f(u,tₙ+h))
    return fzero(g, yₙ + h*f(yₙ,tₙ))
end

#Tests:
function test_linear(integrator::ODEIntegrator)
    @testset "Testing against a linear solution" begin
        #Test 1
        f(y, t) = 1.0
        y₀ = 1.0
        T = 1.0
        h = 0.1
        problem = ODEProblem(f, y₀, T, h)
        solution = solve(problem, integrator)
        @test solution.y ≈ y₀ .+ f(solution.t, 0.0) .* solution.t atol = 1e-10
    end
end

end