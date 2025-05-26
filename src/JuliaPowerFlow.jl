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
export create_ieee9_system, case9, run_power_flow
export run_power_flow_with_visualization

# Add new exports
export convert_to_legacy_format, solve_power_system
export run_power_flow_new, build_ybus_new

# Include core functionality from organized directories
include("../data/case9.jl")
include("../algorithms/power_flow.jl")
include("../visualization/power_flow_visualization.jl")

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

"""
    convert_to_legacy_format(sys::PowerSystem)

Convert PowerSystem to legacy MATPOWER-style data format for compatibility with existing algorithms.
"""
function convert_to_legacy_format(sys::PowerSystem)
    # Initialize legacy format dictionary
    case_data = Dict{String, Any}()
    case_data["version"] = "2"
    case_data["baseMVA"] = sys.base_mva
    
    # Convert buses
    n_buses = length(sys.buses)
    bus_data = zeros(n_buses, 13)
    
    bus_ids = sort(collect(keys(sys.buses)))
    for (i, bus_id) in enumerate(bus_ids)
        bus = sys.buses[bus_id]
        bus_data[i, 1] = bus_id                    # Bus number
        bus_data[i, 2] = Int(bus.bus_type)         # Bus type
        bus_data[i, 3] = bus.load_p                # Pd (MW)
        bus_data[i, 4] = bus.load_q                # Qd (MVAr)
        bus_data[i, 5] = bus.shunt_g               # Gs (MW)
        bus_data[i, 6] = bus.shunt_b               # Bs (MVAr)
        bus_data[i, 7] = 1                         # Area
        bus_data[i, 8] = bus.voltage_magnitude     # Vm (p.u.)
        bus_data[i, 9] = rad2deg(bus.voltage_angle) # Va (degrees)
        bus_data[i, 10] = bus.base_voltage         # Base kV
        bus_data[i, 11] = 1                        # Zone
        bus_data[i, 12] = bus.voltage_max          # Vmax (p.u.)
        bus_data[i, 13] = bus.voltage_min          # Vmin (p.u.)
    end
    case_data["bus"] = bus_data
    
    # Convert generators
    n_gens = length(sys.generators)
    if n_gens > 0
        gen_data = zeros(n_gens, 21)
        gen_ids = sort(collect(keys(sys.generators)))
        
        for (i, gen_id) in enumerate(gen_ids)
            gen = sys.generators[gen_id]
            gen_data[i, 1] = gen.bus_id              # Bus number
            gen_data[i, 2] = gen.p_output            # Pg (MW)
            gen_data[i, 3] = gen.q_output            # Qg (MVAr)
            gen_data[i, 4] = gen.q_max               # Qmax (MVAr)
            gen_data[i, 5] = gen.q_min               # Qmin (MVAr)
            gen_data[i, 6] = gen.voltage_setpoint    # Vg (p.u.)
            gen_data[i, 7] = sys.base_mva           # mBase (MVA)
            gen_data[i, 8] = gen.is_online ? 1 : 0   # Status
            gen_data[i, 9] = gen.p_max               # Pmax (MW)
            gen_data[i, 10] = gen.p_min              # Pmin (MW)
        end
        case_data["gen"] = gen_data
    else
        case_data["gen"] = zeros(0, 21)
    end
    
    # Convert branches
    n_branches = length(sys.branches)
    if n_branches > 0
        branch_data = zeros(n_branches, 13)
        branch_ids = sort(collect(keys(sys.branches)))
        
        for (i, branch_id) in enumerate(branch_ids)
            branch = sys.branches[branch_id]
            branch_data[i, 1] = branch.from_bus       # From bus
            branch_data[i, 2] = branch.to_bus         # To bus
            branch_data[i, 3] = branch.resistance     # R (p.u.)
            branch_data[i, 4] = branch.reactance      # X (p.u.)
            branch_data[i, 5] = branch.susceptance    # B (p.u.)
            branch_data[i, 6] = branch.rate_a         # Rate A (MVA)
            branch_data[i, 7] = branch.rate_b         # Rate B (MVA)
            branch_data[i, 8] = branch.rate_c         # Rate C (MVA)
            branch_data[i, 9] = branch.tap_ratio      # Tap ratio
            branch_data[i, 10] = branch.phase_shift   # Phase shift (degrees)
            branch_data[i, 11] = branch.is_online ? 1 : 0  # Status
        end
        case_data["branch"] = branch_data
    else
        case_data["branch"] = zeros(0, 13)
    end
    
    return case_data
end

"""
    solve_power_system(sys::PowerSystem; max_iter=30, tol=1e-6, verbose=true)

Solve power flow for PowerSystem using existing algorithms.
"""
function solve_power_system(sys::PowerSystem; max_iter=30, tol=1e-6, verbose=true)
    # Convert to legacy format
    case_data = convert_to_legacy_format(sys)
    
    # Use existing power flow solver
    try
        V, S, Sf, St = run_power_flow(case_data)
        
        # Convert results back to PowerSystem format
        results = Dict{String, Any}()
        results["voltage"] = V
        results["power_injection"] = S
        results["line_flows_from"] = Sf
        results["line_flows_to"] = St
        results["converged"] = true
        results["case_data"] = case_data
        
        return results
    catch e
        if verbose
            println("Power flow solution failed: $e")
        end
        return Dict("converged" => false, "error" => string(e))
    end
end

"""
    run_power_flow_new(sys::PowerSystem; with_visualization=false, kwargs...)

Run power flow analysis on PowerSystem with optional visualization.
"""
function run_power_flow_new(sys::PowerSystem; with_visualization=false, kwargs...)
    case_data = convert_to_legacy_format(sys)
    
    if with_visualization
        try
            # Use visualization version if available
            return run_power_flow_with_visualization(case_data)
        catch e
            println("Visualization not available, using standard solver: $e")
            return run_power_flow(case_data)
        end
    else
        return run_power_flow(case_data)
    end
end

"""
    build_ybus_new(sys::PowerSystem)

Build admittance matrix for PowerSystem using existing algorithm.
"""
function build_ybus_new(sys::PowerSystem)
    case_data = convert_to_legacy_format(sys)
    return build_ybus(case_data)
end

"""
    update_system_from_results!(sys::PowerSystem, results::Dict)

Update PowerSystem with power flow results.
"""
function update_system_from_results!(sys::PowerSystem, results::Dict)
    if !results["converged"]
        @warn "Cannot update system: power flow did not converge"
        return
    end
    
    V = results["voltage"]
    case_data = results["case_data"]
    
    # Update bus voltages
    bus_ids = sort(collect(keys(sys.buses)))
    for (i, bus_id) in enumerate(bus_ids)
        if i <= length(V)
            sys.buses[bus_id].voltage_magnitude = abs(V[i])
            sys.buses[bus_id].voltage_angle = angle(V[i])
        end
    end
    
    # Update generator reactive power output from results
    if haskey(results, "case_data") && haskey(case_data, "gen")
        gen_data = case_data["gen"]
        gen_ids = sort(collect(keys(sys.generators)))
        
        # Calculate generator reactive power from power flow results
        for (i, gen_id) in enumerate(gen_ids)
            if i <= size(gen_data, 1)
                # This would require recalculating from the voltage solution
                # For now, we'll leave Q output as specified
            end
        end
    end
end

"""
    analyze_system_performance(sys::PowerSystem)

Analyze the performance of all generators in the system.
"""
function analyze_system_performance(sys::PowerSystem)
    println("="^60)
    println("POWER SYSTEM PERFORMANCE ANALYSIS")
    println("="^60)
    
    # Solve power flow
    results = solve_power_system(sys, verbose=false)
    
    if !results["converged"]
        println("❌ Power flow did not converge!")
        return
    end
    
    println("✅ Power flow converged successfully")
    
    # Update system with results
    update_system_from_results!(sys, results)
    
    # Analyze each generator
    println("\nGenerator Analysis:")
    println("-"^40)
    
    for (gen_id, gen) in sys.generators
        if !gen.is_online
            continue
        end
        
        bus = sys.buses[gen.bus_id]
        println("Generator $gen_id (Bus $(gen.bus_id)):")
        @printf("  Voltage: %.4f ∠ %.2f° p.u.\n", 
                bus.voltage_magnitude, rad2deg(bus.voltage_angle))
        @printf("  Power Output: %.2f MW + j%.2f MVAr\n", 
                gen.p_output, gen.q_output)
        
        # Check capability constraints
        valid, msg = check_capability_constraints(gen, gen.p_output, gen.q_output, bus.voltage_magnitude)
        status = valid ? "✅" : "⚠️"
        println("  Capability Check: $status $msg")
        
        # Calculate Q_min for current P
        q_min_calc = calculate_qmin_function(gen, gen.p_output, bus.voltage_magnitude)
        @printf("  Q_min at current P: %.2f MVAr\n", q_min_calc)
        
        # Calculate utilization
        p_util = gen.p_output / gen.p_max * 100
        s_current = sqrt(gen.p_output^2 + gen.q_output^2)
        s_util = s_current / gen.s_max * 100
        @printf("  Utilization: P=%.1f%%, S=%.1f%%\n", p_util, s_util)
        println()
    end
    
    # System summary
    total_generation_p = sum(gen.p_output for gen in values(sys.generators) if gen.is_online)
    total_load_p = sum(bus.load_p for bus in values(sys.buses))
    total_generation_q = sum(gen.q_output for gen in values(sys.generators) if gen.is_online)
    total_load_q = sum(bus.load_q for bus in values(sys.buses))
    
    println("System Summary:")
    println("-"^40)
    @printf("Total Generation: %.2f MW + j%.2f MVAr\n", total_generation_p, total_generation_q)
    @printf("Total Load: %.2f MW + j%.2f MVAr\n", total_load_p, total_load_q)
    @printf("Active Power Balance: %.4f MW\n", total_generation_p - total_load_p)
    @printf("Reactive Power Balance: %.4f MVAr\n", total_generation_q - total_load_q)
    
    return results
end

end # module JuliaPowerFlow
