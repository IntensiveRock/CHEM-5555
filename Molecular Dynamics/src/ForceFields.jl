module ForceFields

using ..CellLists
export ForceField, compute_pair_interaction, compute_forces!, update_energy!

abstract type ForceField end

function compute_forces!(ff::ForceField)
    r = ff.base.atoms.r
    f = ff.base.atoms.f
    N = ff.base.N
    L = ff.base.L
    
    V = W = 0.0
    f .= 0.0
    rij = zeros(3)
    rcut² = isnothing(ff.parameters.rcut) ? Inf : ff.parameters.rcut^2
    if isa(ff.cells, CellList)
        update_cells!(ff.cells, r, L)
        for i in 1:N
            ri = r[:,i]
            for j in get_cell_neighbors(ff.cells, i, r, L, rcut²)
                rij .= ri - r[:,j]
                rij .= rij .- L .* round.(Int64, rij ./ L)
                r2 = sum(rij .* rij)
                if r2 < rcut²
                    fij, vij = compute_pair_interaction(ff, r2)
                    f[:,i] .+= fij.*rij
                    f[:,j] .-= fij.*rij
                    V += vij
                    W₀ = sum(fij.*rij)
                    W += i < j ? W₀ : -W₀
                end
            end
        end
    else
        for i in 1:(N-1)
            ri = r[:,i]
            for j in (i+1):N
                rij .= ri - r[:,j]
                rij .= rij .- L .* round.(Int64, rij ./ L)
                r2 = sum(rij .* rij)
                if r2 < rcut²
                    fij, vij = compute_pair_interaction(ff, r2)
                    f[:,i] .+= fij.*rij
                    f[:,j] .-= fij.*rij
                    V += vij
                    W += sum(fij.*rij)
                end
            end
        end
    end
    ff.potential_energy = V
    ff.virial = 2.0*W
    return
end

# Interface for force field implementations
function compute_pair_interaction(ff::ForceField, r2::Float64)
    error("compute_pair_interaction not implemented for $(typeof(ff))")
end

function update_energy!(ff::ForceField)
    error("update_energy! not implemented for $(typeof(ff))")
end

end