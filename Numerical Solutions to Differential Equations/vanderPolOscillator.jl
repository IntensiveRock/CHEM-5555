module vanderPolOscillator

include("JDE_ODE_Base.jl")
include("JDE_ODE.jl")
import .JDE_ODE_Base: ODEIntegrator, ODEProblem, step
import .JDE_ODE: Heun, RK4

using Test, Plots

#Define the van der Pol oscillator problem using the structures (struct) and functions defined in JDE_ODE.jl.


end
