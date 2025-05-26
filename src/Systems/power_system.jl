"""
    PowerSystem

Main container for all power system components.
"""
mutable struct PowerSystem
    base_mva::Float64
    buses::Dict{Int, Bus}
    branches::Dict{Int, Branch}
    generators::Dict{Int, Generator}
    loads::Dict{Int, Load}
    
    function PowerSystem(base_mva::Float64 = 100.0)
        new(base_mva, Dict{Int, Bus}(), Dict{Int, Branch}(), 
            Dict{Int, Generator}(), Dict{Int, Load}())
    end
end

"""
    add_component!(sys::PowerSystem, component)

Add a component to the power system.
"""
function add_component!(sys::PowerSystem, bus::Bus)
    sys.buses[bus.id] = bus
end

function add_component!(sys::PowerSystem, branch::Branch)
    sys.branches[branch.id] = branch
end

function add_component!(sys::PowerSystem, generator::Generator)
    sys.generators[generator.id] = generator
end

function add_component!(sys::PowerSystem, load::Load)
    sys.loads[load.id] = load
end

"""
    get_bus_count(sys::PowerSystem)

Get the number of buses in the system.
"""
function get_bus_count(sys::PowerSystem)
    return length(sys.buses)
end

"""
    identify_bus_types(sys::PowerSystem)

Identify PV, PQ, and reference buses.
"""
function identify_bus_types(sys::PowerSystem)
    pv = Int[]
    pq = Int[]
    ref = Int[]
    
    for (id, bus) in sys.buses
        if bus.bus_type == SLACK_BUS
            push!(ref, id)
        elseif bus.bus_type == PV_BUS
            push!(pv, id)
        elseif bus.bus_type == PQ_BUS
            push!(pq, id)
        end
    end
    
    return ref, pv, pq
end
