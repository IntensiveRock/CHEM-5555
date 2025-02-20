module ForceFields

abstract type ForceField end

# Generic force calculation loop that works for any force field
function compute_forces!(sys, ff::ForceField)
    r = sys.base.atoms.r
    f = sys.base.atoms.f
    N = sys.base.N
    L = sys.base.L
    
    V = W = 0.0
    f .= 0.0
    
    for i ∈ 1:N
        ri = r[:,i]
        for j ∈ i+1:N
            rij = ri - r[:, j]
            rij .-= L .* round.(rij ./ L)
            r2 = dot(rij, rij)
            if r2 < ff.rcut^2
                fmag, vij = compute_pair_interaction(ff, r2)  # This is what changes per force field
                fij = fmag*rij
                f[:, i] .+= fij
                f[:, j] .-= fij
                V += vij
                W += dot(fij, rij)
            end
        end
    end
    
    sys.potential_energy = V
    sys.virial = W
end

# Interface that needs to be implemented for each force field
function compute_pair_interaction(ff::ForceField, r2::Float64)
    error("compute_pair_interaction not implemented for $(typeof(ff))")
end

function update_energy!(ff::ForceField)
    error("update_energy! not implemented for $(typeof(ff))")
end

export ForceField, compute_forces!
end