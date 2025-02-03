module ODE_Integrators

import UnPack: @unpack

include("JDE_ODE_Base.jl")
using .JDE_ODE_Base
import .JDE_ODE_Base: step

export Heun, RK2, RK4

struct Heun <: ODEIntegrator end
struct RK2 <: ODEIntegrator end
function step(problem::ODEProblem, ::Union{Heun,RK2}, yₙ::typename, tₙ::Float64) where typename
    @unpack f, h = problem
    k₁ = f(yₙ, tₙ)
    k₂ = f(yₙ + h*k₁, tₙ + h)
    return yₙ + (h/2)*(k₁ + k₂)
end

struct RK4 <: ODEIntegrator end
function step(problem::ODEProblem, ::RK4, yₙ::typename, tₙ::Float64) where typename
    @unpack f, h = problem
    k₁ = f(yₙ, tₙ)
    k₂ = f(yₙ + (h/2)*k₁, tₙ + h/2)
    k₃ = f(yₙ + (h/2)*k₂, tₙ + h/2)
    k₄ = f(yₙ + h*k₃, tₙ + h)
    return yₙ + (h/6)*(k₁ + 2*k₂ + 2*k₃ + k₄)
end

end