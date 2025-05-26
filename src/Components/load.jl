"""
    Load

Represents a load in the power system.
"""
mutable struct Load
    id::Int
    bus_id::Int
    
    # Power consumption
    p_load::Float64        # MW
    q_load::Float64        # MVAr
    
    # Load model parameters
    voltage_dependency_p::Float64  # voltage exponent for P
    voltage_dependency_q::Float64  # voltage exponent for Q
    
    # Status
    is_online::Bool
    
    function Load(id::Int, bus_id::Int;
                 p_load::Float64 = 0.0,
                 q_load::Float64 = 0.0,
                 voltage_dependency_p::Float64 = 0.0,  # constant power
                 voltage_dependency_q::Float64 = 0.0,  # constant power
                 is_online::Bool = true)
        new(id, bus_id, p_load, q_load, voltage_dependency_p, voltage_dependency_q, is_online)
    end
end

"""
    calculate_load_power(load::Load, voltage::Float64)

Calculate the actual load power considering voltage dependency.
"""
function calculate_load_power(load::Load, voltage::Float64)
    if !load.is_online
        return 0.0, 0.0
    end
    
    p_actual = load.p_load * (voltage ^ load.voltage_dependency_p)
    q_actual = load.q_load * (voltage ^ load.voltage_dependency_q)
    
    return p_actual, q_actual
end
