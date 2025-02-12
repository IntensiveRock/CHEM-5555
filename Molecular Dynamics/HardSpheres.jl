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


function packing_fraction_to_density(packing_fraction::Float64)
    return 6*packing_fraction/π
end

function packing_fraction(density::Float64)
    return π*density/6
end

mutable struct hard_spheres
    system::MonoAtomic.system
    dr::Float64
    acceptance::Float64
    function hard_spheres(sys::MonoAtomic.system)
        sys_ = deepcopy(sys)
        MonoAtomic.initalize!(sys_)
        dr_ = 0.01
        acceptance_ = 0.5
        new(sys_, dr_, acceptance_)
    end
    function hard_spheres(filename::String)
        sys_ = MonoAtomic.system(filename)
        if sys_.T == 0.0
            ρ_ = packing_fraction_to_density(sys_.ρ)
            sys_ = MonoAtomic.system(N=sys_.N, ρ=ρ_, T=0.0)
        end
        acceptance_ = 0.5
        MonoAtomic.initalize!(sys_)
        dr_ = 0.01
        new(sys_, dr_, acceptance_)
    end

    function hard_spheres(; N::Int64=1000, ϕ::Float64=0.8)
        ρ_ = packing_fraction_to_density(ϕ)
        sys_ = MonoAtomic.system(N=N, ρ =ρ_, T = 0.0)
        MonoAtomic.initalize!(sys_)
        dr_ = 0.01
        acceptance_ = 0.5
        new(sys_, dr_, acceptance_)
    end
end

end