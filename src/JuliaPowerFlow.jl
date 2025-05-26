module JuliaPowerFlow

using LinearAlgebra, SparseArrays, Printf
using Plots, ColorSchemes, Statistics

# Export main types and functions
export PowerSystem, Bus, Branch, Generator, Load
export BusType, GeneratorType, SLACK_BUS, PV_BUS, PQ_BUS
export THERMAL, HYDRO, WIND, SOLAR, NUCLEAR
export build_ybus, newton_raphson_power_flow
export add_component!, get_bus_count, identify_bus_types
export check_capability_constraints, calculate_qmin_function
export create_ieee9_system, case9

# Include existing functionality from the standalone files
include("../case9.jl")
include("../power_flow.jl")

# Define the component types directly in this module for now
@enum BusType begin
    PQ_BUS = 1
    PV_BUS = 2
    SLACK_BUS = 3
end

@enum GeneratorType begin
    THERMAL
    HYDRO
    WIND
    SOLAR
    NUCLEAR
end

# Bus component definition
mutable struct Bus
    id::Int
    bus_type::BusType
    voltage_magnitude::Float64
    voltage_angle::Float64
    load_p::Float64
    load_q::Float64
    shunt_g::Float64
    shunt_b::Float64
    base_voltage::Float64
    voltage_max::Float64
    voltage_min::Float64
    
    function Bus(id::Int; 
                 bus_type::BusType = PQ_BUS,
                 voltage_magnitude::Float64 = 1.0,
                 voltage_angle::Float64 = 0.0,
                 load_p::Float64 = 0.0,
                 load_q::Float64 = 0.0,
                 shunt_g::Float64 = 0.0,
                 shunt_b::Float64 = 0.0,
                 base_voltage::Float64 = 100.0,
                 voltage_max::Float64 = 1.1,
                 voltage_min::Float64 = 0.9)
        new(id, bus_type, voltage_magnitude, voltage_angle, 
            load_p, load_q, shunt_g, shunt_b, base_voltage, voltage_max, voltage_min)
    end
end

# Generator component definition
mutable struct Generator
    id::Int
    bus_id::Int
    generator_type::GeneratorType
    p_output::Float64
    q_output::Float64
    voltage_setpoint::Float64
    p_max::Float64
    p_min::Float64
    q_max::Float64
    q_min::Float64
    s_max::Float64
    eq_min::Float64
    eq_max::Float64
    xd::Float64
    delta_max::Float64
    is_online::Bool
    
    function Generator(id::Int, bus_id::Int;
                      generator_type::GeneratorType = THERMAL,
                      p_output::Float64 = 0.0,
                      q_output::Float64 = 0.0,
                      voltage_setpoint::Float64 = 1.0,
                      p_max::Float64 = 100.0,
                      p_min::Float64 = 0.0,
                      q_max::Float64 = 50.0,
                      q_min::Float64 = -50.0,
                      s_max::Float64 = 120.0,
                      eq_min::Float64 = 0.8,
                      eq_max::Float64 = 1.4,
                      xd::Float64 = 0.8,
                      delta_max::Float64 = 60.0,
                      is_online::Bool = true)
        new(id, bus_id, generator_type, p_output, q_output, voltage_setpoint,
            p_max, p_min, q_max, q_min, s_max, eq_min, eq_max, xd, delta_max, is_online)
    end
end

# Branch component definition
mutable struct Branch
    id::Int
    from_bus::Int
    to_bus::Int
    resistance::Float64
    reactance::Float64
    susceptance::Float64
    tap_ratio::Float64
    phase_shift::Float64
    rate_a::Float64
    rate_b::Float64
    rate_c::Float64
    is_online::Bool
    
    function Branch(id::Int, from_bus::Int, to_bus::Int;
                   resistance::Float64 = 0.0,
                   reactance::Float64 = 0.1,
                   susceptance::Float64 = 0.0,
                   tap_ratio::Float64 = 1.0,
                   phase_shift::Float64 = 0.0,
                   rate_a::Float64 = 999.0,
                   rate_b::Float64 = 999.0,
                   rate_c::Float64 = 999.0,
                   is_online::Bool = true)
        new(id, from_bus, to_bus, resistance, reactance, susceptance,
            tap_ratio, phase_shift, rate_a, rate_b, rate_c, is_online)
    end
end

# Load component definition
mutable struct Load
    id::Int
    bus_id::Int
    p_load::Float64
    q_load::Float64
    voltage_dependency_p::Float64
    voltage_dependency_q::Float64
    is_online::Bool
    
    function Load(id::Int, bus_id::Int;
                 p_load::Float64 = 0.0,
                 q_load::Float64 = 0.0,
                 voltage_dependency_p::Float64 = 0.0,
                 voltage_dependency_q::Float64 = 0.0,
                 is_online::Bool = true)
        new(id, bus_id, p_load, q_load, voltage_dependency_p, voltage_dependency_q, is_online)
    end
end

# PowerSystem definition
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

# Component management functions
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

function get_bus_count(sys::PowerSystem)
    return length(sys.buses)
end

# Generator capability functions
function check_capability_constraints(gen::Generator, p::Float64, q::Float64, ut::Float64 = 1.0)
    # Apparent power constraint
    if p^2 + q^2 > gen.s_max^2
        return false, "Apparent power limit exceeded"
    end
    
    # Active power limits
    if p < gen.p_min || p > gen.p_max
        return false, "Active power limit exceeded"
    end
    
    # Reactive power limits
    if q < gen.q_min || q > gen.q_max
        return false, "Reactive power limit exceeded"
    end
    
    # Check delta constraints if applicable
    if p > 0
        delta_rad = gen.delta_max * π / 180
        eq_for_p = p * gen.xd / (ut * sin(delta_rad))
        
        if eq_for_p < gen.eq_min || eq_for_p > gen.eq_max
            return false, "Internal voltage or power angle limit exceeded"
        end
    end
    
    return true, "Within capability region"
end

function calculate_qmin_function(gen::Generator, p::Float64, ut::Float64 = 1.0)
    if p <= 0
        return gen.q_min
    end
    
    delta_rad = gen.delta_max * π / 180
    eq_min = gen.eq_min
    
    q_from_delta = (eq_min * ut * cos(delta_rad) - ut^2) / gen.xd
    
    return max(gen.q_min, q_from_delta)
end

# IEEE 9-bus system creation
function create_ieee9_system()
    sys = PowerSystem(100.0)
    
    # Add buses
    add_component!(sys, Bus(1, bus_type=SLACK_BUS, voltage_magnitude=1.04, base_voltage=16.5))
    add_component!(sys, Bus(2, bus_type=PV_BUS, voltage_magnitude=1.025, base_voltage=18.0))
    add_component!(sys, Bus(3, bus_type=PV_BUS, voltage_magnitude=1.025, base_voltage=13.8))
    add_component!(sys, Bus(4, bus_type=PQ_BUS, base_voltage=230.0))
    add_component!(sys, Bus(5, bus_type=PQ_BUS, load_p=125.0, load_q=50.0, base_voltage=230.0))
    add_component!(sys, Bus(6, bus_type=PQ_BUS, load_p=90.0, load_q=30.0, base_voltage=230.0))
    add_component!(sys, Bus(7, bus_type=PQ_BUS, base_voltage=230.0))
    add_component!(sys, Bus(8, bus_type=PQ_BUS, load_p=100.0, load_q=35.0, base_voltage=230.0))
    add_component!(sys, Bus(9, bus_type=PQ_BUS, base_voltage=230.0))
    
    # Add generators
    add_component!(sys, Generator(1, 1, p_output=71.64, voltage_setpoint=1.04, 
                                 p_max=250.0, p_min=10.0, s_max=300.0))
    add_component!(sys, Generator(2, 2, p_output=163.0, voltage_setpoint=1.025,
                                 p_max=300.0, p_min=10.0, s_max=350.0))
    add_component!(sys, Generator(3, 3, p_output=85.0, voltage_setpoint=1.025,
                                 p_max=270.0, p_min=10.0, s_max=320.0))
    
    # Add branches
    add_component!(sys, Branch(1, 1, 4, resistance=0.0, reactance=0.0576, susceptance=0.0))
    add_component!(sys, Branch(2, 4, 5, resistance=0.017, reactance=0.092, susceptance=0.158))
    add_component!(sys, Branch(3, 5, 6, resistance=0.039, reactance=0.17, susceptance=0.358))
    add_component!(sys, Branch(4, 3, 6, resistance=0.0, reactance=0.0586, susceptance=0.0))
    add_component!(sys, Branch(5, 6, 7, resistance=0.0119, reactance=0.1008, susceptance=0.209))
    add_component!(sys, Branch(6, 7, 8, resistance=0.0085, reactance=0.072, susceptance=0.149))
    add_component!(sys, Branch(7, 8, 2, resistance=0.0, reactance=0.0625, susceptance=0.0))
    add_component!(sys, Branch(8, 8, 9, resistance=0.032, reactance=0.161, susceptance=0.306))
    add_component!(sys, Branch(9, 9, 4, resistance=0.01, reactance=0.085, susceptance=0.176))
    
    return sys
end

end # module
