module RDF

#Compute the radial distribution function. And keep a running average.

include("HardSpheres.jl")
using .HardSpheres

function rdf(hs::hard_spheres)
    #Compute the radial distribution function.
    r = positions(hs)
    @show r
end

end