"""
    Generator

Represents a generator in the power system with capability constraints.
"""
mutable struct Generator
    id::Int
    bus_id::Int
    generator_type::GeneratorType
    
    # Power output
    p_output::Float64      # MW
    q_output::Float64      # MVAr
    voltage_setpoint::Float64  # p.u.
    
    # Capability limits
    p_max::Float64         # MW
    p_min::Float64         # MW
    q_max::Float64         # MVAr
    q_min::Float64         # MVAr
    s_max::Float64         # MVA
    
    # Synchronous machine parameters
    eq_min::Float64        # Internal voltage min (p.u.)
    eq_max::Float64        # Internal voltage max (p.u.)
    xd::Float64           # Direct-axis reactance (p.u.)
    delta_max::Float64    # Maximum power angle (degrees)
    
    # Status
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

"""
    check_capability_constraints(gen::Generator, p::Float64, q::Float64, ut::Float64 = 1.0)

Check if the given P and Q are within the generator capability region.
"""
function check_capability_constraints(gen::Generator, p::Float64, q::Float64, ut::Float64 = 1.0)
    # Apparent power constraint
    if p^2 + q^2 > gen.s_max^2
        return false, "Apparent power limit exceeded"
    end
    
    # Active power limits
    if p < gen.p_min || p > gen.p_max
        return false, "Active power limit exceeded"
    end
    
    # Reactive power limits (constant)
    if q < gen.q_min || q > gen.q_max
        return false, "Reactive power limit exceeded"
    end
    
    # Check delta constraints if applicable
    if p > 0  # Only check for positive power
        delta_rad = gen.delta_max * π / 180
        eq_for_p = p * gen.xd / (ut * sin(delta_rad))
        
        if eq_for_p < gen.eq_min || eq_for_p > gen.eq_max
            return false, "Internal voltage or power angle limit exceeded"
        end
    end
    
    return true, "Within capability region"
end

"""
    calculate_qmin_function(gen::Generator, p::Float64, ut::Float64 = 1.0)

Calculate Q_min as a function of P for the generator.
"""
function calculate_qmin_function(gen::Generator, p::Float64, ut::Float64 = 1.0)
    if p <= 0
        return gen.q_min
    end
    
    delta_rad = gen.delta_max * π / 180
    eq_min = gen.eq_min
    
    # Q = (Eq * Ut * cos(delta) - Ut^2) / Xd
    q_from_delta = (eq_min * ut * cos(delta_rad) - ut^2) / gen.xd
    
    return max(gen.q_min, q_from_delta)
end
