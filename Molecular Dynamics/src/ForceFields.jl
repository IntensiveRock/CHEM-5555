module ForceFields

abstract type ForceField end

export ForceField, compute_pair_interaction, compute_forces!, update_energy!
# Generic force calculation loop that works for any force field
function compute_forces!(ff::ForceField)
    r = ff.base.atoms.r
    f = ff.base.atoms.f
    N = ff.base.N
    L = ff.base.L
    
    V = W = 0.0
    f .= 0.0
    @inline dot(x,y) = sum(x .* y)
    rij = zeros(3)
    ri = zeros(3)
    rcut² = isnothing(ff.parameters.rcut) ? Inf : ff.parameters.rcut^2
    for i ∈ 1:N
        ri = r[:,i]
        for j ∈ (i+1):N
            rij .= ri - r[:, j]
            rij .= rij .- L .* round.(Int64, rij ./ L)
            r2 = dot(rij, rij)
            if r2 < rcut²
                fmag, vij = compute_pair_interaction(ff, r2)  # This is what changes per force field
                fij = fmag*rij
                f[:, i] .+= fij
                f[:, j] .-= fij
                V += vij
                W += dot(fij, rij)
            end
        end
    end
    ff.potential_energy = V
    ff.virial = 2.0*W
end

# Interface that needs to be implemented for each force field
function compute_pair_interaction(ff::ForceField, r2::Float64)
    error("compute_pair_interaction not implemented for $(typeof(ff))")
end

function update_energy!(ff::ForceField)
    error("update_energy! not implemented for $(typeof(ff))")
end

end