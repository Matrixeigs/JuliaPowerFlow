"""
    Bus

Represents a bus in the power system.
"""
mutable struct Bus
    id::Int
    bus_type::BusType
    voltage_magnitude::Float64  # p.u.
    voltage_angle::Float64      # radians
    load_p::Float64            # MW
    load_q::Float64            # MVAr
    shunt_g::Float64           # MW at 1 p.u. voltage
    shunt_b::Float64           # MVAr at 1 p.u. voltage
    base_voltage::Float64      # kV
    voltage_max::Float64       # p.u.
    voltage_min::Float64       # p.u.
    
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

"""
    get_voltage_complex(bus::Bus)

Get the complex voltage of the bus.
"""
function get_voltage_complex(bus::Bus)
    return bus.voltage_magnitude * exp(im * bus.voltage_angle)
end
