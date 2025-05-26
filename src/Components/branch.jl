"""
    Branch

Represents a transmission line or transformer in the power system.
"""
mutable struct Branch
    id::Int
    from_bus::Int
    to_bus::Int
    
    # Electrical parameters
    resistance::Float64     # p.u.
    reactance::Float64      # p.u.
    susceptance::Float64    # p.u.
    
    # Transformer parameters
    tap_ratio::Float64      # transformer tap ratio
    phase_shift::Float64    # phase shift angle (degrees)
    
    # Ratings
    rate_a::Float64         # MVA rating
    rate_b::Float64         # MVA rating (emergency)
    rate_c::Float64         # MVA rating (emergency)
    
    # Status
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

"""
    get_admittance(branch::Branch)

Calculate the series admittance of the branch.
"""
function get_admittance(branch::Branch)
    z = complex(branch.resistance, branch.reactance)
    return 1.0 / z
end

"""
    get_shunt_admittance(branch::Branch)

Calculate the shunt admittance of the branch.
"""
function get_shunt_admittance(branch::Branch)
    return complex(0.0, branch.susceptance)
end
