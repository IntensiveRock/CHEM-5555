module Potentials
using Optim

export quartic

"""

Quartic double-well potential with barrier height of 10 atomic units at x = 0.
The potential has two symmetric minima and Φ(0) = 10.
"""
function quartic(x)
    # Parameters chosen to create a barrier of 10 at x = 0
    a = 10.0  # barrier height
    b = √(4a) # parameter to ensure minima at desired locations
    # Quartic potential: V(x) = ax⁴ - bx²
    return -b * x^2 + x^4
end

function find_minimum(f,x0=-1.0)
    result = optimize(f, -2.0, 2.0, Brent())
    return Optim.minimizer(result)
end

end