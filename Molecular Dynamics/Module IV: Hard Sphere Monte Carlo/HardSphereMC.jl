module HardSphereMC

include("../src/MonoAtomic.jl")
include("HardSpheres.jl")
include("../src/Accessors.jl")
include("../src/Calculations.jl")

include("../src/MSD.jl")
include("../src/RDF.jl")

include("../src/MonteCarlo.jl")

using .MonoAtomic
using .HardSpheres
using .Accessors
using .Calculations

using .MSD
using .RDF

using .MonteCarlo

#Must include the structs that you use for the calculations.
export msd, rdf
# Export the types and functions needed
export hard_spheres, sweep!, mc_simulation!, box_length, positions, number_of_particles
export initialize!, compute!, finalize!

end