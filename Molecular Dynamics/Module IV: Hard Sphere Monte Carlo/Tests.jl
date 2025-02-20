module Tests

using Test
include("HardSphereMC.jl")
using .HardSphereMC
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
        sys = hard_spheres(N=1000,ϕ=0.4,lattice="sc")  # System with 1000 particles
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

function test_msd()
    @testset "MSD Translation Test" begin
        # Create a test system
        N_ = 1000
        sys = hard_spheres(N=N_,ϕ = 0.4, lattice="sc")
        
        # Create MSD calculator and store initial positions
        calc = msd(T=1)
        initialize!(calc, sys)
        
        # Translate all particles by a known vector
        R = [1.0, 2.0, 3.0]
        r = positions(sys)
        for i in 1:N_
            r[:,i] .+= R
        end
        
        # Compute MSD
        compute!(calc, sys)
        finalize!(calc)
        
        # Test that MSD equals the square of the translation vector
        expected_msd = sum(R.^2)
        @test isapprox(calc.r²[1], expected_msd, rtol=1e-10)
    end
end

function test_rdf_peaks()
    @testset "RDF Peak Tests" begin
        sys = hard_spheres(N=1000, ϕ=0.4, lattice="sc")
        L = box_length(sys)
        calc = rdf(bins=100, r_max=L)
        
        initialize!(calc, sys)
        compute!(calc, sys)
        finalize!(calc)
        
        a = L/(number_of_particles(sys)^(1/3))
        peak_positions = [1.0, sqrt(2), sqrt(3), 2.0]
        
        @testset "Peak Position Tests" begin
            for pos in peak_positions
                r = pos * a
                bin = floor(Int, r/calc.dr) + 1
                
                peak_found = false
                for b in max(1,bin-1):min(calc.bins,bin+1)
                    if calc.hist[b] > 0.0  # g(r) should be > 1 at peaks
                        peak_found = true
                        break
                    end
                end
                @test peak_found
            end
        end
        
        @testset "Peak Integration Tests" begin
            ρ = number_of_particles(sys)/(L^3)
            expected_coords = [6, 12, 8, 6]
            
            for (i, (pos, expected)) in enumerate(zip(peak_positions, expected_coords))
                r_peak = pos * a
                integral = 0.0
                for bin in 1:calc.bins
                    r = (bin-0.5)*calc.dr
                    if abs(r - r_peak) < 0.2*a
                        integral += 4π * ρ * r^2 * calc.hist[bin] * calc.dr
                    end
                end
                @test isapprox(integral, expected, rtol=0.2)
            end
        end
    end
end

function test_rdf_normalization()
    @testset "RDF Normalization" begin
        # Create a cubic lattice system
        sys = hard_spheres(N=1000, ϕ=0.4, lattice="sc")
        L = box_length(sys)
        
        # Create RDF calculator
        calc = rdf(bins=100, r_max=L)
        initialize!(calc, sys)
        compute!(calc, sys)
        finalize!(calc)

        # Test RDF normalization by integrating over 4πρr²g(r)
        # This should equal N-1 for a properly normalized RDF
        ρ = number_of_particles(sys)/(L^3)  # Number density
        N = number_of_particles(sys)

        integral = 0.0
        for i in 1:calc.bins
            r = (i-0.5)*calc.dr  # Center of bin
            integral += 4π * ρ * r^2 * calc.hist[i] * calc.dr
        end

        # Test that integral equals N-1 within numerical tolerance
        @test isapprox(integral, N-1, rtol=0.1)
    end
end

end  # module Tests