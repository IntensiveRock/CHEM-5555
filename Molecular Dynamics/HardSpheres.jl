module HardSpheres
#A Monte Carlo simulation of the hard sphere fluid

include("MonoAtomic.jl")
using .MonoAtomic

export hard_spheres
# Usage here would be the following:
# sys = MonoAtomic.system(N=1000, ϕ=0.5)
# hs = hard_spheres(sys)

#Or:
# hs = hard_spheres("parameters.toml")

#Or:
# hs = hard_spheres(N=1000, ϕ=0.5)

mutable struct hard_spheres
    sys::MonoAtomic.system
    dr::Float64
    acceptance::Float64
    r::Array{Float64,2}
    function hard_spheres(sys::MonoAtomic.system)
        sys_ = sys
        initalize!(sys_)
        dr_ = 0.1
        acceptance_ = 0.0
        r_ = @views sys.r #The analogy of passing by reference in Julia
        new(sys, dr_, acceptance_, r_)
    end
    function hard_spheres(filename::String)
        sys_ = MonoAtomic.system(filename)
        initalize!(sys_)
        dr_ = 0.1
        acceptance_ = 0.0
        r_ = @views sys.r
        new(sys_, dr_, acceptance_, r_)
    end
    function hard_spheres(; N::Int64=1000, ϕ::Float64=0.8)
        sys_ = MonoAtomic.system(N=N, ϕ=ϕ)
        initalize!(sys_)
        dr_ = 0.1
        acceptance_ = 0.0
        r_ = @views sys.r
        new(sys_, dr_, acceptance_, r_)
    end
end

end