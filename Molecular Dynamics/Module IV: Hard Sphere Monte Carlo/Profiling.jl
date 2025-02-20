module Profiling

using StatProfilerHTML
using BenchmarkTools
include("HardSphereMC.jl")
using .HardSphereMC

export profile, time_it

function profile(; N=343, ϕ=0.4,steps=1000)
    hs = hard_spheres(N=N, ϕ=ϕ)
    @time sweep!(hs)
    @profilehtml begin 
        for i in 1:steps
            sweep!(hs)
        end
    end
end

function time_it(; N=343, ϕ=0.4,steps=1000)
    hs = hard_spheres(N=N, ϕ=ϕ)
    @btime sweep!($hs) 
end

end