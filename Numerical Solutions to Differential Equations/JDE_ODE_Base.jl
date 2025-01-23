module JDE_ODE_Base

export ODEIntegrator, ODEProblem, ODESolution, solve, step

import UnPack: @unpack

abstract type ODEIntegrator end

struct ODESolution{typename}
    t::Vector{Float64}
    y::Vector{typename}
end

struct ODEProblem{typename}
    f::Function #Forcing function f = f(y,t)
    y₀::typename #Initial condition
    T::Float64 #Final time
    h::Float64 #Time step
end

#This function must be overwritten in the integrator implementation.
function step end

#Implementation of the solve function. It is a fairly general function that can be used with any ODEProblem and ODEIntegrator.
function solve(problem::ODEProblem{typename},integrator::ODEIntegrator) where typename
    @unpack h, T, y₀ = problem
    t = 0.0:h:T
    N = length(t)
    y = Vector{typename}(undef,N)
    y[1] = y₀
    for n ∈ 1:(N-1)
        y[n+1] = step(problem, integrator, y[n], t[n])
    end
    return ODESolution{typename}( collect(t), y)
end

end