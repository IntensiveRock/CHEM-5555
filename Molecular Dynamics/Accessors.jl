module Accessors

using Test
include("HardSpheres.jl")
using .HardSpheres

export positions, velocities, forces, box_length, number_of_particles

function get_field(system::typename, x::Symbol) where typename
    function search_nested(obj)
        # If we find the field directly, return it
        if hasproperty(obj, x)
            return getfield(obj, x)
        end
        
        # Otherwise, search through all fields recursively
        for field in fieldnames(typeof(obj))
            if isdefined(obj, field) && !isprimitivetype(typeof(getfield(obj, field)))
                field_value = getfield(obj, field)
                try
                    return search_nested(field_value)
                catch
                    continue
                end
            end
        end
        throw(ErrorException("Field not found"))
    end

    try
        return search_nested(system)
    catch
        fields = fieldnames(typename)
        error("Could not find field '$x' in system of type $(typename) or its members. Available fields: $(fields)")
    end
end

function positions(system::typename) where typename
    get_field(system, :r)
end

function velocities(system::typename) where typename
    get_field(system, :v)
end

function forces(system::typename) where typename
    get_field(system, :f)
end

function box_length(system::typename) where typename
    get_field(system, :L)
end

function number_of_particles(system::typename) where typename
    get_field(system, :N)
end


end