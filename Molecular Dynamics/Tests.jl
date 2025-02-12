module Tests

using Test
include("HardSpheres.jl")
using .HardSpheres
include("Accessors.jl")
using .Accessors

#Test accessor functions.
function test_positions()
    @testset "Accessors Tests" begin
        # Create a test system
        sys = hard_spheres()  # System with 1000 particles
    
        @test positions(sys) === sys.system.atoms.r #Ensure that the return is equivalent at the level of the memory address.
    end
end

function test_lattice()
    @testset "Lattice Tests" begin
        # Create a test system
        sys = hard_spheres()  # System with 1000 particles
        L = box_length(sys)
        N = number_of_particles(sys)
        r = positions(sys)

        a = L/N^(1/3)  # Lattice constant
        expected_positions = zeros(Float64, 3, N)
        # Initialize the positions of the particles in a cubic lattice using a different, more obscure method.
        for i in 1:N
            x = a * ((i-1) % N^(1/3))
            y = a * (((i-1) ÷ N^(1/3)) % N^(1/3))
            z = a * (((i-1) ÷ (N^(1/3)^2)) % N^(1/3))
            expected_positions[:, i] .= [x, y, z]
        end
    
        # Check that the positions are initialized correctly
        expected_positions = sys.system.atoms.r
        @test r == expected_positions 

        # Check that the positions are within the box length
        @testset "Testing box boundaries" begin
            for i in 1:N
                @test r[1, i] >= 0 && r[1, i] <= L
                @test r[2, i] >= 0 && r[2, i] <= L
                @test r[3, i] >= 0 && r[3, i] <= L
            end
        end
        # Check that the positions are unique
        @testset "Testing uniqueness of positions" begin
            for i in 1:N
                for j in i+1:N
                    @test r[:, i] != r[:, j]
                end
            end
        end
        # Check to make sure that a is the minimum separation distance between particles.
        @testset "Testing minimum separation distance" begin
            a² = Inf  # Initialize with Inf to find minimum
            for i in 1:N
                for j in (i+1):N
                    rij = r[:, i] - r[:, j]
                    rij = rij .- L * round.(rij ./ L)  # Apply periodic boundary conditions without Int64 conversion
                    d² = sum(rij .* rij)
                    if d² < 1e-10  # Debug zero distances
                        println("Zero distance found between particles $i and $j")
                        println("r[$i] = ", r[:, i])
                        println("r[$j] = ", r[:, j])
                        println("rij = ", rij)
                    end
                    a² = min(a², d²)  # Keep track of minimum d²
                end
            end
            @test isapprox(sqrt(a²), a)
        end
    end
end

end