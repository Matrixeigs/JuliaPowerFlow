using LinearAlgebra, DifferentialEquations, Plots, SparseArrays

"""
Electromagnetic Transients simulation for fault current analysis in AC distribution systems
with Inverter-Based Distributed Energy Resources (IBDERs)
"""

struct ExternalGridParameters
    voltage_magnitude::Float64    # Grid voltage magnitude (V)
    frequency::Float64           # Grid frequency (Hz)
    impedance::ComplexF64        # Grid impedance (Ohms)
    inertia::Float64            # Grid inertia constant (s)
    damping::Float64            # Damping coefficient
    time_constant::Float64       # Voltage regulation time constant (s)
end

struct IBDERParameters
    rated_power::Float64      # Rated power (VA)
    rated_voltage::Float64    # Rated voltage (V)
    max_current::Float64      # Maximum current (A)
    current_limit::Float64    # Current limit factor (pu)
    response_time::Float64    # Control response time (s)
    voltage_support::Bool     # Voltage support capability
    frequency_support::Bool   # Frequency support capability
    # New parameters for ODE/DAE model
    filter_inductance::Float64   # Output filter inductance (H)
    filter_capacitance::Float64  # Output filter capacitance (F)
    dc_voltage::Float64         # DC bus voltage (V)
    current_controller_kp::Float64  # Current controller proportional gain
    current_controller_ki::Float64  # Current controller integral gain
end

struct FaultParameters
    fault_type::Symbol       # :three_phase, :single_line_ground, :line_to_line
    fault_resistance::Float64 # Fault resistance (Ohms)
    fault_location::Int      # Bus number where fault occurs
    fault_time::Float64      # Fault initiation time (s)
    fault_duration::Float64  # Fault duration (s)
end

struct NetworkData
    bus_count::Int
    line_impedances::Matrix{ComplexF64}  # Line impedance matrix
    bus_voltages::Vector{ComplexF64}     # Pre-fault bus voltages
    ibder_locations::Vector{Int}         # Buses with IBDERs
    ibder_params::Vector{IBDERParameters}
    load_impedances::Vector{ComplexF64}  # Load impedances at each bus
    grid_bus::Int                        # External grid connection bus
    grid_params::ExternalGridParameters  # External grid parameters
end

# State vector indices for ODE system
struct StateIndices
    grid_angle::Int
    grid_frequency::Int
    grid_voltage_ref::Int
    ibder_current_d::Vector{Int}
    ibder_current_q::Vector{Int}
    ibder_voltage_d::Vector{Int}
    ibder_voltage_q::Vector{Int}
    ibder_integrator_d::Vector{Int}
    ibder_integrator_q::Vector{Int}
end

function simulate_emt_fault_current(network::NetworkData, fault::FaultParameters, 
                                  simulation_time::Float64, dt::Float64=1e-5)
    """
    Main EMT simulation function with ODE/DAE models for grid and IBDERs
    """
    
    # Initialize state vector
    n_ibders = length(network.ibder_locations)
    n_states = 3 + 6*n_ibders
    
    # Create state indices structure
    state_idx = create_state_indices(n_ibders)
    
    # Initial conditions with proper steady-state values
    u0 = zeros(n_states)
    u0[state_idx.grid_angle] = 0.0
    u0[state_idx.grid_frequency] = 2π * network.grid_params.frequency
    u0[state_idx.grid_voltage_ref] = network.grid_params.voltage_magnitude
    
    # Initialize IBDER states with proper steady-state values
    for i in 1:n_ibders
        # Initial current based on rated power (positive for generation)
        i_rated = network.ibder_params[i].rated_power / (sqrt(3) * network.ibder_params[i].rated_voltage)
        u0[state_idx.ibder_current_d[i]] = i_rated  # Positive d-axis current for real power
        u0[state_idx.ibder_current_q[i]] = 0.0      # Zero q-axis current initially
        u0[state_idx.ibder_voltage_d[i]] = network.ibder_params[i].rated_voltage
        u0[state_idx.ibder_voltage_q[i]] = 0.0
        u0[state_idx.ibder_integrator_d[i]] = 0.0
        u0[state_idx.ibder_integrator_q[i]] = 0.0
    end
    
    # Build admittance matrix
    Y_bus = build_admittance_matrix(network)
    
    # Define ODE system with autodiff-compatible operations
    function ode_system!(du, u, p, t)
        network, fault, Y_bus, state_idx = p
        
        # Initialize du to zero
        fill!(du, 0.0)
        
        # Check if fault is active
        fault_active = (t >= fault.fault_time && 
                       t <= fault.fault_time + fault.fault_duration)
        
        # Grid dynamics
        grid_ode!(du, u, network.grid_params, state_idx, fault_active)
        
        # IBDER dynamics
        for i in 1:length(network.ibder_locations)
            ibder_ode!(du, u, network.ibder_params[i], state_idx, i, t, fault_active)
        end
        
        # Network algebraic equations (simplified for autodiff compatibility)
        solve_network_algebraic_autodiff!(du, u, network, Y_bus, state_idx, fault, fault_active)
    end
    
    # Solve ODE system with finite differences to avoid autodiff issues
    tspan = (0.0, simulation_time)
    prob = ODEProblem(ode_system!, u0, tspan, (network, fault, Y_bus, state_idx))
    
    println("Starting EMT simulation with finite difference solver...")
    # Use FBDF solver which is less sensitive to Jacobian accuracy
    sol = solve(prob, FBDF(autodiff=false), dt=dt, adaptive=true, reltol=1e-4, abstol=1e-6)
    
    # Extract results
    return extract_emt_results(sol, network, fault, state_idx)
end

function create_state_indices(n_ibders::Int)
    """Create state indices for the ODE system"""
    idx = 1
    
    # Grid states
    grid_angle = idx; idx += 1
    grid_frequency = idx; idx += 1
    grid_voltage_ref = idx; idx += 1
    
    # IBDER states
    ibder_current_d = Int[]
    ibder_current_q = Int[]
    ibder_voltage_d = Int[]
    ibder_voltage_q = Int[]
    ibder_integrator_d = Int[]
    ibder_integrator_q = Int[]
    
    for i in 1:n_ibders
        push!(ibder_current_d, idx); idx += 1
        push!(ibder_current_q, idx); idx += 1
        push!(ibder_voltage_d, idx); idx += 1
        push!(ibder_voltage_q, idx); idx += 1
        push!(ibder_integrator_d, idx); idx += 1
        push!(ibder_integrator_q, idx); idx += 1
    end
    
    return StateIndices(grid_angle, grid_frequency, grid_voltage_ref,
                       ibder_current_d, ibder_current_q, ibder_voltage_d,
                       ibder_voltage_q, ibder_integrator_d, ibder_integrator_q)
end

function grid_ode!(du, u, grid_params::ExternalGridParameters, state_idx::StateIndices, fault_active::Bool)
    """External grid ODE model with proper power feedback"""
    
    # Extract grid states
    θ = u[state_idx.grid_angle]
    ω = u[state_idx.grid_frequency]
    V_ref = u[state_idx.grid_voltage_ref]
    
    # Grid dynamics
    du[state_idx.grid_angle] = ω - 2π * grid_params.frequency
    
    # Frequency dynamics with power feedback
    ω_base = 2π * grid_params.frequency
    Δω = ω - ω_base
    
    # Simplified power calculation (would be more accurate with network solution)
    P_mech = 1.0  # Base mechanical power to maintain frequency
    
    # During fault, electrical power drops, causing frequency rise
    # After fault, system should return to nominal frequency
    if fault_active
        P_elec = 0.5 * P_mech  # Reduced electrical power during fault
    else
        P_elec = P_mech  # Normal electrical power
    end
    
    # Governor/frequency response
    P_ref_adjustment = -0.05 * Δω  # Droop response
    P_net = P_mech + P_ref_adjustment - P_elec
    
    # Swing equation: dω/dt = (P_net - D*Δω) / (2H)
    du[state_idx.grid_frequency] = (P_net - grid_params.damping * Δω) / (2 * grid_params.inertia)
    
    # Voltage reference dynamics (AVR response)
    V_set = grid_params.voltage_magnitude
    
    if fault_active
        # During fault, voltage reference tries to maintain voltage
        du[state_idx.grid_voltage_ref] = (V_set * 1.1 - V_ref) / grid_params.time_constant
    else
        # After fault, return to nominal voltage
        du[state_idx.grid_voltage_ref] = (V_set - V_ref) / grid_params.time_constant
    end
end

function solve_network_algebraic_autodiff!(du, u, network::NetworkData, Y_bus::Matrix{ComplexF64}, 
                                         state_idx::StateIndices, fault::FaultParameters, fault_active::Bool)
    """Solve network algebraic equations with proper voltage calculation"""
    
    # Modify admittance matrix for fault
    Y_fault_real = zeros(2*network.bus_count, 2*network.bus_count)
    Y_fault = fault_active ? apply_fault_admittance(Y_bus, fault) : Y_bus
    
    # Convert complex admittance matrix to real equivalent form
    for i in 1:network.bus_count
        for j in 1:network.bus_count
            Y_fault_real[i, j] = real(Y_fault[i, j])
            Y_fault_real[i, j + network.bus_count] = -imag(Y_fault[i, j])
            Y_fault_real[i + network.bus_count, j] = imag(Y_fault[i, j])
            Y_fault_real[i + network.bus_count, j + network.bus_count] = real(Y_fault[i, j])
        end
    end
    
    # Current injections in real form
    I_injection_real = zeros(2*network.bus_count)
    
    # Grid current injection - improved grid model
    θ_grid = u[state_idx.grid_angle]
    V_grid = u[state_idx.grid_voltage_ref]
    
    grid_bus = network.grid_bus
    if grid_bus <= network.bus_count
        # Grid voltage phasor
        V_grid_real = V_grid * cos(θ_grid)
        V_grid_imag = V_grid * sin(θ_grid)
        
        # Grid behind impedance model
        Z_grid = network.grid_params.impedance
        Y_grid = 1.0 / Z_grid
        
        # Add grid admittance to the matrix
        Y_fault_real[grid_bus, grid_bus] += real(Y_grid)
        Y_fault_real[grid_bus, grid_bus + network.bus_count] += -imag(Y_grid)
        Y_fault_real[grid_bus + network.bus_count, grid_bus] += imag(Y_grid)
        Y_fault_real[grid_bus + network.bus_count, grid_bus + network.bus_count] += real(Y_grid)
        
        # Grid voltage source injection
        I_grid_real = real(Y_grid) * V_grid_real - imag(Y_grid) * V_grid_imag
        I_grid_imag = imag(Y_grid) * V_grid_real + real(Y_grid) * V_grid_imag
        
        I_injection_real[grid_bus] += I_grid_real
        I_injection_real[grid_bus + network.bus_count] += I_grid_imag
    end
    
    # IBDER current injections - corrected signs
    for (i, bus_idx) in enumerate(network.ibder_locations)
        if bus_idx <= network.bus_count
            i_d = u[state_idx.ibder_current_d[i]]
            i_q = u[state_idx.ibder_current_q[i]]
            
            # Correct dq to real/imaginary conversion for current injection
            # Positive i_d should inject real current, positive i_q should inject reactive current
            I_injection_real[bus_idx] += i_d
            I_injection_real[bus_idx + network.bus_count] += i_q  # Corrected sign
        end
    end
    
    # Solve for bus voltages
    try
        reg_factor = 1e-8
        Y_reg = Y_fault_real + reg_factor * I(2*network.bus_count)
        V_solution = Y_reg \ I_injection_real
        
        # Update IBDER terminal voltage states
        for (i, bus_idx) in enumerate(network.ibder_locations)
            if bus_idx <= network.bus_count
                v_real = V_solution[bus_idx]
                v_imag = V_solution[bus_idx + network.bus_count]
                
                # Fast voltage tracking
                time_constant = 0.001
                du[state_idx.ibder_voltage_d[i]] = (v_real - u[state_idx.ibder_voltage_d[i]]) / time_constant
                du[state_idx.ibder_voltage_q[i]] = (v_imag - u[state_idx.ibder_voltage_q[i]]) / time_constant
            end
        end
    catch e
        println("Warning: Network solution failed, using previous values")
        for (i, bus_idx) in enumerate(network.ibder_locations)
            du[state_idx.ibder_voltage_d[i]] = 0.0
            du[state_idx.ibder_voltage_q[i]] = 0.0
        end
    end
end

function ibder_ode!(du, u, ibder_param::IBDERParameters, state_idx::StateIndices, 
                   ibder_idx::Int, t::Float64, fault_active::Bool)
    """IBDER ODE/DAE model with corrected current signs and references"""
    
    # Extract IBDER states
    i_d = u[state_idx.ibder_current_d[ibder_idx]]
    i_q = u[state_idx.ibder_current_q[ibder_idx]]
    v_d = u[state_idx.ibder_voltage_d[ibder_idx]]
    v_q = u[state_idx.ibder_voltage_q[ibder_idx]]
    xi_d = u[state_idx.ibder_integrator_d[ibder_idx]]
    xi_q = u[state_idx.ibder_integrator_q[ibder_idx]]
    
    # System parameters
    L = ibder_param.filter_inductance
    V_dc = ibder_param.dc_voltage
    kp = ibder_param.current_controller_kp
    ki = ibder_param.current_controller_ki
    ω = 2π * 60.0
    
    # Current references based on operating mode
    if fault_active
        # During fault: current limiting and voltage support
        i_max = ibder_param.max_current * ibder_param.current_limit
        
        if ibder_param.voltage_support
            # Voltage magnitude for voltage support logic
            v_mag = sqrt(v_d^2 + v_q^2)
            if v_mag > 0.1 * ibder_param.rated_voltage
                # Reduce active current, increase reactive current for voltage support
                i_d_ref = min(0.5 * i_max, ibder_param.rated_power / (2 * sqrt(3) * ibder_param.rated_voltage))
                # Positive q-axis current for reactive power injection (voltage support)
                i_q_available = sqrt(max(0.0, i_max^2 - i_d_ref^2))
                i_q_ref = min(i_q_available, 0.8 * i_max)
            else
                i_d_ref = 0.0
                i_q_ref = 0.0
            end
        else
            # Current limiting only - maintain positive active current
            i_d_ref = min(i_max, ibder_param.rated_power / (sqrt(3) * ibder_param.rated_voltage))
            i_q_ref = 0.0
        end
    else
        # Normal operation - positive d-axis current for real power generation
        i_d_ref = ibder_param.rated_power / (sqrt(3) * ibder_param.rated_voltage)
        i_q_ref = 0.0
    end
    
    # Current controller errors
    e_d = i_d_ref - i_d
    e_q = i_q_ref - i_q
    
    # PI controller for current control with proper decoupling
    v_d_ref = kp * e_d + ki * xi_d - ω * L * i_q
    v_q_ref = kp * e_q + ki * xi_q + ω * L * i_d
    
    # PWM modulation with DC voltage limit
    v_dc_half = V_dc / 2
    m_d = max(-0.95, min(0.95, v_d_ref / v_dc_half))
    m_q = max(-0.95, min(0.95, v_q_ref / v_dc_half))
    
    # Inverter output voltages
    v_inv_d = m_d * v_dc_half
    v_inv_q = m_q * v_dc_half
    
    # Current dynamics (inductor equation) - corrected signs
    R = 0.01
    du[state_idx.ibder_current_d[ibder_idx]] = (v_inv_d - v_d - R*i_d + ω*L*i_q) / L
    du[state_idx.ibder_current_q[ibder_idx]] = (v_inv_q - v_q - R*i_q - ω*L*i_d) / L
    
    # Integrator states for PI controllers
    du[state_idx.ibder_integrator_d[ibder_idx]] = e_d
    du[state_idx.ibder_integrator_q[ibder_idx]] = e_q
end

function extract_emt_results(sol, network::NetworkData, fault::FaultParameters, state_idx::StateIndices)
    """Extract results from ODE solution with proper transformations"""
    
    t = sol.t
    n_steps = length(t)
    n_ibders = length(network.ibder_locations)
    
    # Initialize result arrays
    bus_voltages = zeros(ComplexF64, network.bus_count, n_steps)
    ibder_currents = zeros(ComplexF64, n_ibders, n_steps)
    grid_frequency = zeros(Float64, n_steps)
    grid_voltage = zeros(Float64, n_steps)
    
    for i in 1:n_steps
        u = sol.u[i]
        
        # Extract grid states
        grid_frequency[i] = u[state_idx.grid_frequency] / (2π)
        grid_voltage[i] = u[state_idx.grid_voltage_ref]
        
        # Extract IBDER currents with correct transformation
        for j in 1:n_ibders
            i_d = u[state_idx.ibder_current_d[j]]
            i_q = u[state_idx.ibder_current_q[j]]
            # Correct dq to phasor transformation
            ibder_currents[j, i] = i_d + 1im * i_q  # Corrected sign
            
            # Reconstruct bus voltages
            bus_idx = network.ibder_locations[j]
            if bus_idx <= network.bus_count
                v_d = u[state_idx.ibder_voltage_d[j]]
                v_q = u[state_idx.ibder_voltage_q[j]]
                bus_voltages[bus_idx, i] = v_d + 1im * v_q  # Corrected transformation
            end
        end
        
        # Add grid bus voltage
        if network.grid_bus <= network.bus_count
            θ_grid = u[state_idx.grid_angle]
            V_grid = u[state_idx.grid_voltage_ref]
            bus_voltages[network.grid_bus, i] = V_grid * exp(1im * θ_grid)
        end
    end
    
    return EMTResultsODE(t, bus_voltages, ibder_currents, grid_frequency, grid_voltage, fault)
end

struct EMTResultsODE
    time::Vector{Float64}
    bus_voltages::Matrix{ComplexF64}
    ibder_currents::Matrix{ComplexF64}
    grid_frequency::Vector{Float64}
    grid_voltage::Vector{Float64}
    fault_params::FaultParameters
end

function build_admittance_matrix(network::NetworkData)
    """Build the network admittance matrix"""
    Y = zeros(ComplexF64, network.bus_count, network.bus_count)
    
    # Add line admittances
    for i in 1:network.bus_count
        for j in 1:network.bus_count
            if i != j && network.line_impedances[i, j] != 0
                y_ij = 1.0 / network.line_impedances[i, j]
                Y[i, j] = -y_ij
                Y[i, i] += y_ij
            end
        end
    end
    
    # Add load admittances
    for i in 1:network.bus_count
        if network.load_impedances[i] != 0
            Y[i, i] += 1.0 / network.load_impedances[i]
        end
    end
    
    return Y
end

function apply_fault_admittance(Y_bus::Matrix{ComplexF64}, fault::FaultParameters)
    """Apply fault admittance to the system"""
    Y_fault = copy(Y_bus)
    fault_bus = fault.fault_location
    
    if fault.fault_type == :three_phase
        # Three-phase fault: add fault admittance to ground
        y_fault = 1.0 / (fault.fault_resistance + 1e-6)  # Small resistance to avoid singularity
        Y_fault[fault_bus, fault_bus] += y_fault
    elseif fault.fault_type == :single_line_ground
        # Single line-to-ground fault
        y_fault = 1.0 / (fault.fault_resistance + 1e-6)
        Y_fault[fault_bus, fault_bus] += y_fault
    elseif fault.fault_type == :line_to_line
        # Line-to-line fault (simplified representation)
        y_fault = 1.0 / (fault.fault_resistance + 1e-6)
        Y_fault[fault_bus, fault_bus] += y_fault * 0.866  # √3/2 factor
    end
    
    return Y_fault
end

function initialize_ibder_states(ibder_params::Vector{IBDERParameters})
    """Initialize IBDER control states"""
    states = []
    for param in ibder_params
        state = Dict(
            :current_ref => 0.0 + 0.0im,
            :voltage_ref => param.rated_voltage,
            :active_power => param.rated_power,
            :reactive_power => 0.0,
            :current_limit_active => false,
            :fault_detected => false
        )
        push!(states, state)
    end
    return states
end

function update_ibder_response(network::NetworkData, states::Vector, 
                             bus_voltages::Vector{ComplexF64}, time::Float64, 
                             dt::Float64, fault_active::Bool)
    """Update IBDER control response during simulation"""
    currents = zeros(ComplexF64, length(network.ibder_locations))
    
    for (idx, bus_idx) in enumerate(network.ibder_locations)
        param = network.ibder_params[idx]
        state = states[idx]
        
        # Get bus voltage
        v_bus = bus_voltages[bus_idx]
        v_mag = abs(v_bus)
        
        # Fault detection logic
        voltage_threshold = 0.7 * param.rated_voltage  # 70% voltage dip threshold
        state[:fault_detected] = v_mag < voltage_threshold || fault_active
        
        if state[:fault_detected]
            # Current limiting during fault
            i_max = param.max_current * param.current_limit
            
            # Voltage support if enabled
            if param.voltage_support && v_mag > 0.1 * param.rated_voltage
                # Provide reactive current for voltage support
                q_ref = min(0.5 * param.rated_power, 
                           i_max * v_mag * 0.5)  # 50% reactive current
                state[:reactive_power] = q_ref
            else
                state[:reactive_power] = 0.0
            end
            
            # Calculate current considering limits
            s_apparent = sqrt(state[:active_power]^2 + state[:reactive_power]^2)
            if s_apparent > 0
                current_mag = min(i_max, s_apparent / v_mag)
                phase_angle = angle(v_bus) - acos(state[:active_power] / s_apparent)
                currents[idx] = current_mag * exp(1im * phase_angle)
            end
            
            state[:current_limit_active] = abs(currents[idx]) >= i_max * 0.95
        else
            # Normal operation
            state[:active_power] = param.rated_power
            state[:reactive_power] = 0.0
            
            if v_mag > 0.1 * param.rated_voltage
                i_magnitude = param.rated_power / v_mag
                currents[idx] = i_magnitude * exp(1im * (angle(v_bus) - π/2))
            end
        end
    end
    
    return currents, states
end

function calculate_branch_currents(network::NetworkData, voltages::Vector{ComplexF64}, 
                                 Y_bus::Matrix{ComplexF64})
    """Calculate branch currents from bus voltages"""
    currents = zeros(ComplexF64, network.bus_count, network.bus_count)
    
    for i in 1:network.bus_count
        for j in 1:network.bus_count
            if i != j && network.line_impedances[i, j] != 0
                # Current from bus i to bus j
                currents[i, j] = (voltages[i] - voltages[j]) / network.line_impedances[i, j]
            end
        end
    end
    
    return currents
end

function plot_emt_results(results::EMTResultsODE; bus_to_plot::Int=1)
    """Plot EMT simulation results for ODE-based simulation"""
    
    # Plot bus voltage magnitude
    p1 = plot(results.time * 1000, abs.(results.bus_voltages[bus_to_plot, :]),
              title="Bus Voltage Magnitude", xlabel="Time (ms)", ylabel="Voltage (V)",
              linewidth=2, legend=false)
    
    # Plot IBDER current magnitude
    if size(results.ibder_currents, 1) > 0
        p2 = plot(results.time * 1000, abs.(results.ibder_currents[1, :]),
                  title="IBDER Current Magnitude", xlabel="Time (ms)", ylabel="Current (A)",
                  linewidth=2, legend=false, color=:red)
    else
        p2 = plot(title="No IBDER Current Data", xlabel="Time (ms)", ylabel="Current (A)")
    end
    
    # Plot grid frequency
    p3 = plot(results.time * 1000, results.grid_frequency,
              title="Grid Frequency", xlabel="Time (ms)", ylabel="Frequency (Hz)",
              linewidth=2, legend=false, color=:green)
    
    # Plot grid voltage reference
    p4 = plot(results.time * 1000, results.grid_voltage,
              title="Grid Voltage Reference", xlabel="Time (ms)", ylabel="Voltage (V)",
              linewidth=2, legend=false, color=:blue)
    
    # Add fault timing indicators
    fault_start = results.fault_params.fault_time * 1000
    fault_end = (results.fault_params.fault_time + results.fault_params.fault_duration) * 1000
    
    vline!(p1, [fault_start, fault_end], color=:red, linestyle=:dash, alpha=0.5, label="Fault")
    vline!(p2, [fault_start, fault_end], color=:red, linestyle=:dash, alpha=0.5, label="Fault")
    vline!(p3, [fault_start, fault_end], color=:red, linestyle=:dash, alpha=0.5, label="Fault")
    vline!(p4, [fault_start, fault_end], color=:red, linestyle=:dash, alpha=0.5, label="Fault")
    
    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end

# Add a plotting function specifically for the enhanced results
function plot_enhanced_emt_results(results::EMTResultsODE)
    """Enhanced plotting function for ODE-based EMT results"""
    
    # Time in milliseconds
    t_ms = results.time * 1000
    
    # Plot 1: All bus voltage magnitudes
    p1 = plot(title="Bus Voltage Magnitudes", xlabel="Time (ms)", ylabel="Voltage (V)")
    for i in 1:size(results.bus_voltages, 1)
        if any(abs.(results.bus_voltages[i, :]) .> 0)  # Only plot if there's data
            plot!(p1, t_ms, abs.(results.bus_voltages[i, :]), 
                  label="Bus $i", linewidth=2)
        end
    end
    
    # Plot 2: IBDER current components
    p2 = plot(title="IBDER Current Components", xlabel="Time (ms)", ylabel="Current (A)")
    if size(results.ibder_currents, 1) > 0
        plot!(p2, t_ms, real.(results.ibder_currents[1, :]), 
              label="Real Current", linewidth=2, color=:red)
        plot!(p2, t_ms, imag.(results.ibder_currents[1, :]), 
              label="Imaginary Current", linewidth=2, color=:blue)
        plot!(p2, t_ms, abs.(results.ibder_currents[1, :]), 
              label="Magnitude", linewidth=2, color=:black, linestyle=:dash)
    end
    
    # Plot 3: Grid frequency deviation
    p3 = plot(title="Grid Frequency Response", xlabel="Time (ms)", ylabel="Frequency (Hz)")
    nominal_freq = 60.0  # Assuming 60 Hz system
    plot!(p3, t_ms, results.grid_frequency, 
          label="Grid Frequency", linewidth=2, color=:green)
    hline!(p3, [nominal_freq], label="Nominal (60 Hz)", 
           linestyle=:dash, color=:gray, alpha=0.7)
    
    # Plot 4: Grid voltage reference
    p4 = plot(title="Grid Voltage Reference", xlabel="Time (ms)", ylabel="Voltage (V)")
    plot!(p4, t_ms, results.grid_voltage, 
          label="Voltage Reference", linewidth=2, color=:purple)
    
    # Add fault timing indicators to all plots
    fault_start = results.fault_params.fault_time * 1000
    fault_end = (results.fault_params.fault_time + results.fault_params.fault_duration) * 1000
    
    for p in [p1, p2, p3, p4]
        vline!(p, [fault_start], color=:red, linestyle=:dot, alpha=0.8, 
               linewidth=2, label="Fault Start")
        vline!(p, [fault_end], color=:orange, linestyle=:dot, alpha=0.8, 
               linewidth=2, label="Fault End")
    end
    
    return plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 900))
end

# Example usage function
function example_emt_simulation()
    """Example of how to use the enhanced EMT simulation with ODE models"""
    
    # Define IBDER parameters with ODE model parameters
    ibder_param = IBDERParameters(
        100e3, 480.0, 150.0, 1.2, 0.01, true, false,
        1e-3, 10e-6, 800.0, 10.0, 100.0  # L, C, V_dc, kp, ki
    )
    
    # Define external grid parameters
    grid_param = ExternalGridParameters(
        13800.0,    # Voltage magnitude (V)
        60.0,       # Frequency (Hz)  
        0.1+0.5im,  # Grid impedance
        5.0,        # Inertia (s)
        0.1,        # Damping
        0.5         # Time constant (s)
    )
    
    # Define network with external grid
    network = NetworkData(
        3,  # 3-bus system
        [0.0+0.0im 0.1+0.3im 0.0+0.0im;
         0.1+0.3im 0.0+0.0im 0.2+0.4im;
         0.0+0.0im 0.2+0.4im 0.0+0.0im],
        [1.0+0.0im, 0.95-0.1im, 0.9-0.05im],
        [2],  # IBDER at bus 2
        [ibder_param],
        [0.0+0.0im, 1.0+0.5im, 0.8+0.4im],
        1,  # Grid at bus 1
        grid_param
    )
    
    # Define fault
    fault = FaultParameters(:three_phase, 0.01, 1, 0.1, 0.05)
    
    # Run simulation
    println("Running EMT simulation...")
    results = simulate_emt_fault_current(network, fault, 0.5, 1e-5)
    
    # Plot results
    println("Generating plots...")
    p1 = plot_emt_results(results, bus_to_plot=1)
    p2 = plot_enhanced_emt_results(results)
    
    display(p1)
    display(p2)
    
    return results
end

println("Enhanced EMT simulation with ODE/DAE models loaded successfully!")
println("Use example_emt_simulation() to run a demonstration")
