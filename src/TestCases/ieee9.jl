"""
    create_ieee9_system()

Create the IEEE 9-bus test system using the new component structure.
"""
function create_ieee9_system()
    sys = PowerSystem(100.0)  # 100 MVA base
    
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
