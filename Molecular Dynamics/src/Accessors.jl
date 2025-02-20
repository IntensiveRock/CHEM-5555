module Accessors

"""
    @generate_accessor name field
Generates an accessor function that recursively searches for field in an object.
"""
macro generate_accessor(name, field)
    return quote
        function $(esc(name))(system::T) where T
            function find_field(obj)
                if hasproperty(obj, $(QuoteNode(field)))
                    return getfield(obj, $(QuoteNode(field)))
                end
                
                for f in fieldnames(typeof(obj))
                    val = getfield(obj, f)
                    if !isprimitivetype(typeof(val))
                        try 
                            return find_field(val)
                        catch
                            continue
                        end
                    end
                end
                throw(ErrorException("Field $($(QuoteNode(field))) not found"))
            end
            
            try
                find_field(system)
            catch
                error("Could not find field '$($(QuoteNode(field)))' in system of type $T")
            end
        end
    end
end

# Generate all accessors
@generate_accessor positions r
@generate_accessor velocities v
@generate_accessor forces f
@generate_accessor box_length L
@generate_accessor number_of_particles N
@generate_accessor time_step dt


export positions, velocities, forces, box_length, number_of_particles, time_step

end