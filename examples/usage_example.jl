"""
Example Usage of JuliaPowerFlow Module

This file demonstrates the key features and functionality of the JuliaPowerFlow package:
1. Creating power system components
2. Building power systems
3. Running power flow analysis (using existing functions)
4. Analyzing generator capability curves
5. Integration with existing power flow code
"""

# Import the JuliaPowerFlow module
push!(LOAD_PATH, pwd() * "/src")
using JuliaPowerFlow

# Also import required packages for this example
using Printf, Plots

println("="^60)
println("JULIA POWER FLOW MODULE - USAGE EXAMPLE")
println("="^60)

#=
EXAMPLE 1: Creating Individual Components
=#
println("\n1. Creating Power System Components:")
println("-" * "="^40)

# Create buses with different types
bus1 = Bus(1, bus_type=SLACK_BUS, voltage_magnitude=1.04, base_voltage=345.0)
bus2 = Bus(2, bus_type=PV_BUS, voltage_magnitude=1.025, load_p=100.0, load_q=35.0)
bus3 = Bus(3, bus_type=PQ_BUS, load_p=125.0, load_q=50.0)

println("✓ Created 3 buses:")
println("  - Bus 1: Slack bus ($(bus1.voltage_magnitude) p.u.)")
println("  - Bus 2: PV bus with load ($(bus2.load_p) MW)")
println("  - Bus 3: PQ bus with load ($(bus3.load_p) MW)")

# Create generators with capability constraints
gen1 = Generator(1, 1, 
                generator_type=THERMAL,
                p_output=150.0, 
                voltage_setpoint=1.04,
                p_max=250.0, p_min=10.0,
                q_max=100.0, q_min=-50.0,
                s_max=300.0,
                eq_min=0.8, eq_max=1.4,
                xd=0.8, delta_max=60.0)

gen2 = Generator(2, 2,
                generator_type=HYDRO,
                p_output=80.0,
                voltage_setpoint=1.025,
                p_max=150.0, p_min=5.0,
                s_max=180.0)

println("✓ Created 2 generators:")
println("  - Gen 1: Thermal ($(gen1.p_max) MW max)")
println("  - Gen 2: Hydro ($(gen2.p_max) MW max)")

# Create transmission lines
line1 = Branch(1, 1, 2, resistance=0.02, reactance=0.08, susceptance=0.1)
line2 = Branch(2, 2, 3, resistance=0.03, reactance=0.12, susceptance=0.15)

println("✓ Created 2 transmission lines")

#=
EXAMPLE 2: Building a Complete Power System
=#
println("\n2. Building Complete Power System:")
println("-" * "="^40)

# Create power system
sys = PowerSystem(100.0)  # 100 MVA base

# Add all components
add_component!(sys, bus1)
add_component!(sys, bus2)
add_component!(sys, bus3)
add_component!(sys, gen1)
add_component!(sys, gen2)
add_component!(sys, line1)
add_component!(sys, line2)

println("✓ Power system created with:")
println("  - $(get_bus_count(sys)) buses")
println("  - $(length(sys.generators)) generators")
println("  - $(length(sys.branches)) branches")
println("  - $(sys.base_mva) MVA base power")

#=
EXAMPLE 3: Using Existing Power Flow Functions
=#
println("\n3. Power Flow Analysis using Existing Functions:")
println("-" * "="^40)

try
    # Use the existing case9 function and power flow solver
    case_data = case9()
    println("✓ Loaded IEEE 9-bus case from existing case9() function")
    
    # Run power flow using existing functions
    V, S, Sf, St = run_power_flow(case_data)
    
    println("✓ Power flow completed using existing solver")
    
catch e
    println("⚠ Power flow calculation encountered an issue:")
    println("  Error: $e")
end

#=
EXAMPLE 4: Generator Capability Analysis
=#
println("\n4. Generator Capability Analysis:")
println("-" * "="^40)

# Test generator capability constraints
test_points = [
    (50.0, 25.0, "Normal operation"),
    (200.0, 80.0, "High power"),
    (100.0, 120.0, "Overexcited"),
    (250.0, 0.0, "Maximum power"),
    (300.0, 50.0, "Overcapacity")
]

println("✓ Testing generator capability constraints:")
for (p, q, description) in test_points
    valid, msg = check_capability_constraints(gen1, p, q)
    status = valid ? "✓" : "✗"
    @printf("  %s P=%.0f MW, Q=%.0f MVAr (%s): %s\n", 
            status, p, q, description, msg)
end

# Calculate Q_min as function of P
println("\n✓ Q_min function for different power levels:")
for p in [0.0, 50.0, 100.0, 150.0, 200.0]  # Changed to Float64 values
    q_min = calculate_qmin_function(gen1, p)
    @printf("  P = %3.0f MW → Q_min = %.2f MVAr\n", p, q_min)
end

#=
EXAMPLE 5: Using Pre-built Test Cases
=#
println("\n5. IEEE Test Case Example:")
println("-" * "="^40)

# Create IEEE 9-bus system using new structure
ieee9_sys = create_ieee9_system()

println("✓ IEEE 9-bus system created using new structure:")
println("  - $(get_bus_count(ieee9_sys)) buses")
println("  - $(length(ieee9_sys.generators)) generators")
println("  - $(length(ieee9_sys.branches)) branches")

# Show bus types
bus_types = []
for (id, bus) in ieee9_sys.buses
    if bus.bus_type == SLACK_BUS
        push!(bus_types, "Bus $id: Slack")
    elseif bus.bus_type == PV_BUS
        push!(bus_types, "Bus $id: PV")
    else
        push!(bus_types, "Bus $id: PQ")
    end
end

println("✓ Bus classification:")
for bt in bus_types
    println("  - $bt")
end

#=
EXAMPLE 6: Comparison with Existing Implementation
=#
println("\n6. Comparison with Existing Implementation:")
println("-" * "="^40)

# Compare data structures
case_data = case9()
println("✓ Original case9() structure:")
println("  - Buses: $(size(case_data["bus"], 1))")
println("  - Branches: $(size(case_data["branch"], 1))")
println("  - Generators: $(size(case_data["gen"], 1))")

println("✓ New PowerSystem structure:")
println("  - Buses: $(length(ieee9_sys.buses))")
println("  - Branches: $(length(ieee9_sys.branches))")
println("  - Generators: $(length(ieee9_sys.generators))")

#=
EXAMPLE 7: Component Modification
=#
println("\n7. Dynamic Component Modification:")
println("-" * "="^40)

# Modify generator output
original_p = gen1.p_output
gen1.p_output = 180.0
println("✓ Modified Gen 1 output: $(original_p) → $(gen1.p_output) MW")

# Take a component offline
gen2.is_online = false
println("✓ Took Gen 2 offline")

# Modify load
original_load = bus2.load_p
bus2.load_p = 120.0
println("✓ Modified Bus 2 load: $(original_load) → $(bus2.load_p) MW")

# Restore original values
gen1.p_output = original_p
gen2.is_online = true
bus2.load_p = original_load
println("✓ Restored original component settings")

#=
EXAMPLE 8: Integration Test
=#
println("\n8. Integration Test:")
println("-" * "="^40)

# Test using both existing power flow and new components
try
    println("✓ Testing integration between new components and existing solver...")
    
    # Create a simple 3-bus system
    simple_sys = PowerSystem(100.0)
    add_component!(simple_sys, Bus(1, bus_type=SLACK_BUS, voltage_magnitude=1.0))
    add_component!(simple_sys, Bus(2, bus_type=PQ_BUS, load_p=50.0, load_q=20.0))
    add_component!(simple_sys, Generator(1, 1, p_output=60.0, s_max=100.0))
    add_component!(simple_sys, Branch(1, 1, 2, resistance=0.01, reactance=0.1))
    
    println("  - Created simple 2-bus system with new components")
    println("  - System has $(get_bus_count(simple_sys)) buses")
    
    # Test generator capability
    p_test, q_test = 60.0, 25.0
    valid, msg = check_capability_constraints(simple_sys.generators[1], p_test, q_test)
    println("  - Generator capability test: $valid ($msg)")
    
catch e
    println("  Note: Full integration requires converter functions")
    println("  Error: $e")
end

#=
SUMMARY
=#
println("\n" * "="^60)
println("EXAMPLE SUMMARY")
println("="^60)

println("""
The JuliaPowerFlow module provides:

1. COMPONENT MODELING:
   • Bus, Generator, Branch, and Load components
   • Comprehensive parameter specifications
   • Component status management

2. SYSTEM BUILDING:
   • PowerSystem container for all components
   • Easy component addition with add_component!()
   • Automatic system validation

3. COMPATIBILITY:
   • Works with existing case9() and power flow functions
   • Maintains backward compatibility
   • Provides new object-oriented interface

4. GENERATOR CAPABILITY:
   • Full synchronous machine constraints
   • P-Q capability region analysis
   • Dynamic Q_min calculation

5. EXTENSIBILITY:
   • Modular design for easy extension
   • Support for custom components
   • Flexible analysis workflows

To use this module in your own projects:
   using JuliaPowerFlow
   sys = PowerSystem()
   # Add components...
   # Use existing power flow: run_power_flow(case9())
   # Use new components: check_capability_constraints(gen, p, q)
""")

println("Example completed successfully! ✓")