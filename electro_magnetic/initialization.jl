# initialization.jl - System Initialization Functions

# ============================================================================
# COMPONENT INITIALIZATION
# ============================================================================

"""
Initialize all system components with default or specified parameters

Returns:
- components: SystemComponents struct containing all initialized components
"""
function initialize_all_components()
    println("  - Initializing synchronous generator...")
    generator = SyncGenerator()
    
    println("  - Initializing PV system...")
    pv_system = PVSystem()
    
    println("  - Initializing BESS...")
    bess = BESSSystem()
    
    println("  - Initializing transmission line...")
    line = TransmissionLine()
    
    println("  - Initializing transformer...")
    transformer = Transformer()
    
    println("  - Initializing induction motor...")
    motor = InductionMotor()
    
    println("  - Initializing EV charger...")
    ev_charger = EVCharger()
    
    components = SystemComponents(generator, pv_system, bess, line, transformer, motor, ev_charger)
    
    println("  - All components initialized successfully!")
    return components
end

"""
Set initial conditions for all system states

Returns:
- x0: initial state vector
"""
function set_initial_conditions()
    println("  - Setting generator initial conditions...")
    println("  - Setting PV system initial conditions...")
    println("  - Setting BESS initial conditions...")
    println("  - Setting motor initial conditions...")
    println("  - Setting EV charger initial conditions...")
    println("  - Setting network initial conditions...")
    
    x0 = zeros(N_STATES)
    
    # ========================================================================
    # GENERATOR INITIAL CONDITIONS
    # ========================================================================
    
    # Calculate steady-state operating point
    P_gen_init = GEN_PARAMS[:P_gen]  # 0.8 pu
    Q_gen_init = GEN_PARAMS[:Q_gen]  # 0.3 pu
    V_t_init = GEN_PARAMS[:V_t]      # 1.0 pu
    
    # Assume generator connected to infinite bus
    δ_init = asin(P_gen_init / V_t_init)  # Power angle
    
    # DQ currents from power and voltage
    I_gen_mag = sqrt(P_gen_init^2 + Q_gen_init^2) / V_t_init
    φ = atan(Q_gen_init, P_gen_init)
    
    x0[STATE_INDICES[:gen_id]] = I_gen_mag * cos(φ)      # d-axis current
    x0[STATE_INDICES[:gen_iq]] = I_gen_mag * sin(φ)      # q-axis current
    x0[STATE_INDICES[:gen_ifd]] = 1.2                     # Field current
    x0[STATE_INDICES[:gen_ikd]] = 0.0                     # d-axis damper current
    x0[STATE_INDICES[:gen_ikq]] = 0.0                     # q-axis damper current
    x0[STATE_INDICES[:gen_ωr]] = 1.0                      # Rotor speed (pu)
    x0[STATE_INDICES[:gen_δ]] = δ_init                    # Power angle
    
    # ========================================================================
    # PV SYSTEM INITIAL CONDITIONS
    # ========================================================================
    
    # Standard test conditions
    G_init = 1000.0  # W/m²
    T_init = 298.15  # K
    
    # Calculate initial PV operating point
    V_pv_init = PV_PARAMS[:Vmp]
    I_pv_init = PV_PARAMS[:Imp]
    
    # DC-DC converter initial conditions
    d_boost_init = 0.6  # Initial duty cycle
    Vdc_pv_init = V_pv_init / (1 - d_boost_init)
    
    x0[STATE_INDICES[:pv_iL]] = I_pv_init                 # Inductor current
    x0[STATE_INDICES[:pv_Vdc]] = Vdc_pv_init             # DC voltage
    x0[STATE_INDICES[:pv_mppt_int]] = 0.0                # MPPT integrator
    x0[STATE_INDICES[:pv_curr_int_d]] = 0.0              # Current controller d
    x0[STATE_INDICES[:pv_curr_int_q]] = 0.0              # Current controller q
    
    # ========================================================================
    # BESS INITIAL CONDITIONS
    # ========================================================================
    
    # Battery initial conditions
    SOC_init = 0.8  # 80% state of charge
    I_bat_init = 0.0  # Initially no charging/discharging
    
    # Calculate initial battery voltage
    Q_extracted_init = BESS_PARAMS[:Q_max] * (1 - SOC_init)
    V_bat_init = battery_voltage_model(I_bat_init, SOC_init, Q_extracted_init, BESSSystem())
    
    # DC link voltage
    Vdc_bess_init = 800.0  # V
    
    x0[STATE_INDICES[:bess_ibat]] = I_bat_init            # Battery current
    x0[STATE_INDICES[:bess_SOC]] = SOC_init              # State of charge
    x0[STATE_INDICES[:bess_Vdc]] = Vdc_bess_init         # DC voltage
    x0[STATE_INDICES[:bess_curr_int_d]] = 0.0            # Current controller d
    x0[STATE_INDICES[:bess_curr_int_q]] = 0.0            # Current controller q
    x0[STATE_INDICES[:bess_volt_int_d]] = 0.0            # Voltage controller d
    x0[STATE_INDICES[:bess_volt_int_q]] = 0.0            # Voltage controller q
    x0[STATE_INDICES[:bess_ω_int]] = 0.0                 # Frequency integrator
    
    # ========================================================================
    # MOTOR INITIAL CONDITIONS
    # ========================================================================
    
    # Motor steady-state operating point
    # Assume motor running at rated conditions
    slip_init = 0.03  # 3% slip
    ωr_motor_init = (1 - slip_init) * ω_BASE
    
    # Steady-state currents (simplified calculation)
    V_motor_init = MOTOR_PARAMS[:V_rated] / V_BASE  # pu
    
    # Approximate steady-state currents
    I_motor_mag = 0.8  # pu
    power_factor = 0.85
    φ_motor = acos(power_factor)
    
    x0[STATE_INDICES[:motor_ids]] = I_motor_mag * cos(φ_motor)    # Stator d current
    x0[STATE_INDICES[:motor_iqs]] = I_motor_mag * sin(φ_motor)    # Stator q current
    x0[STATE_INDICES[:motor_idr]] = 0.8 * x0[STATE_INDICES[:motor_ids]]  # Rotor d current
    x0[STATE_INDICES[:motor_iqr]] = 0.8 * x0[STATE_INDICES[:motor_iqs]]  # Rotor q current
    x0[STATE_INDICES[:motor_ωr]] = ωr_motor_init             # Rotor speed
    
    # ========================================================================
    # EV CHARGER INITIAL CONDITIONS
    # ========================================================================
    
    # EV battery initial conditions
    SOC_ev_init = 0.5   # 50% charged
    I_charging_init = 0.0  # Initially not charging
    
    x0[STATE_INDICES[:ev_iac]] = I_charging_init          # AC current
    x0[STATE_INDICES[:ev_SOC]] = SOC_ev_init             # State of charge
    x0[STATE_INDICES[:ev_Vdc]] = EV_PARAMS[:E0]          # DC voltage
    x0[STATE_INDICES[:ev_curr_int]] = 0.0                # Current integrator
    x0[STATE_INDICES[:ev_volt_int]] = 0.0                # Voltage integrator
    
    # ========================================================================
    # NETWORK INITIAL CONDITIONS
    # ========================================================================
    
    # Transmission line currents (steady-state power flow)
    I_line_init = 0.8  # pu current
    
    x0[STATE_INDICES[:line_ia]] = I_line_init * cos(0)           # Phase A
    x0[STATE_INDICES[:line_ib]] = I_line_init * cos(-2π/3)      # Phase B
    x0[STATE_INDICES[:line_ic]] = I_line_init * cos(2π/3)       # Phase C
    
    # Bus voltages (balanced three-phase)
    V_bus_init = 1.0  # pu voltage
    
    x0[STATE_INDICES[:bus_va]] = V_bus_init * cos(0)            # Phase A
    x0[STATE_INDICES[:bus_vb]] = V_bus_init * cos(-2π/3)       # Phase B
    x0[STATE_INDICES[:bus_vc]] = V_bus_init * cos(2π/3)        # Phase C
    
    # Transformer flux linkages
    flux_init = V_bus_init / ω_BASE  # Steady-state flux
    
    x0[STATE_INDICES[:trafo_flux_a]] = flux_init * cos(0)       # Phase A flux
    x0[STATE_INDICES[:trafo_flux_b]] = flux_init * cos(-2π/3)  # Phase B flux
    x0[STATE_INDICES[:trafo_flux_c]] = flux_init * cos(2π/3)   # Phase C flux
    
    # Fault current (initially zero)
    x0[STATE_INDICES[:fault_current]] = 0.0
    
    # ========================================================================
    # INITIAL CONDITION VERIFICATION
    # ========================================================================
    
    println("  - Verifying initial conditions...")
    
    # Check power balance
    P_gen_check = x0[STATE_INDICES[:gen_id]] * x0[STATE_INDICES[:bus_va]] + 
                  x0[STATE_INDICES[:gen_iq]] * x0[STATE_INDICES[:bus_vb]]
    P_load_check = 0.8  # Total load
    
    if abs(P_gen_check - P_load_check) > 0.1
        @warn "Power balance check failed: P_gen = $P_gen_check, P_load = $P_load_check"
    end
    
    # Check voltage magnitudes
    V_mag_check = sqrt(x0[STATE_INDICES[:bus_va]]^2 + x0[STATE_INDICES[:bus_vb]]^2 + 
                      x0[STATE_INDICES[:bus_vc]]^2) / √3
    
    if abs(V_mag_check - 1.0) > 0.05
        @warn "Voltage magnitude check failed: V_mag = $V_mag_check"
    end
    
    # Check for any NaN or Inf values
    for i in 1:length(x0)
        if !isfinite(x0[i])
            error("Initial condition x0[$i] is not finite: $(x0[i])")
        end
    end
    
    println("  - Initial conditions set and verified successfully!")
    
    return x0
end

"""
Calculate steady-state operating point for a given component

Parameters:
- component_type: symbol indicating component type
- parameters: component parameters
- operating_conditions: dict of operating conditions

Returns:
- steady_state: dict of steady-state values
"""
function calculate_steady_state(component_type::Symbol, parameters, operating_conditions)
    steady_state = Dict()
    
    if component_type == :generator
        # Synchronous generator steady-state
        P = operating_conditions[:P]
        Q = operating_conditions[:Q]
        V = operating_conditions[:V]
        
        # Power angle
        δ = asin(P / V)
        
        # Terminal current
        I_t = sqrt(P^2 + Q^2) / V
        φ = atan(Q, P)
        
        steady_state[:delta] = δ
        steady_state[:I_d] = I_t * cos(φ)
        steady_state[:I_q] = I_t * sin(φ)
        steady_state[:omega] = 1.0
        
    elseif component_type == :pv
        # PV system steady-state
        G = operating_conditions[:irradiance]
        T = operating_conditions[:temperature]
        
        # Calculate MPP operating point
        V_mpp = parameters[:Vmp] * (1 + parameters[:Kv] * (T - T_REF))
        I_mpp = parameters[:Imp] * (1 + parameters[:Ki] * (T - T_REF)) * G / G_REF
        
        steady_state[:V_pv] = V_mpp
        steady_state[:I_pv] = I_mpp
        steady_state[:P_pv] = V_mpp * I_mpp
        
    elseif component_type == :bess
        # BESS steady-state
        SOC = operating_conditions[:SOC]
        I_bat = operating_conditions[:I_bat]
        
        # initialization.jl - System Initialization Functions (Continued)

        Q_extracted = parameters[:Q_max] * (1 - SOC)
        V_bat = battery_voltage_model(I_bat, SOC, Q_extracted, BESSSystem())
        
        steady_state[:V_bat] = V_bat
        steady_state[:SOC] = SOC
        steady_state[:I_bat] = I_bat
        
    elseif component_type == :motor
        # Induction motor steady-state
        P_mech = operating_conditions[:P_mechanical]
        V_supply = operating_conditions[:V_supply]
        
        # Approximate steady-state slip
        s = 0.03  # 3% slip assumption
        
        # Motor equivalent circuit solution (simplified)
        R_total = parameters[:Rs] + parameters[:Rr] / s
        X_total = parameters[:Xs] + parameters[:Xr]
        Z_total = sqrt(R_total^2 + X_total^2)
        
        I_motor = V_supply / Z_total
        power_factor = R_total / Z_total
        
        steady_state[:slip] = s
        steady_state[:I_motor] = I_motor
        steady_state[:power_factor] = power_factor
        steady_state[:omega_r] = (1 - s) * ω_BASE
        
    end
    
    return steady_state
end

"""
Initialize control system states

Parameters:
- components: system components

Returns:
- control_states: dictionary of control system initial states
"""
function initialize_control_states(components)
    control_states = Dict()
    
    # PLL states
    control_states[:pll_theta] = 0.0
    control_states[:pll_omega] = ω_BASE
    control_states[:pll_integrator] = 0.0
    
    # Current controller states for each inverter
    control_states[:pv_curr_int_d] = 0.0
    control_states[:pv_curr_int_q] = 0.0
    control_states[:bess_curr_int_d] = 0.0
    control_states[:bess_curr_int_q] = 0.0
    control_states[:ev_curr_int] = 0.0
    
    # Voltage controller states
    control_states[:bess_volt_int_d] = 0.0
    control_states[:bess_volt_int_q] = 0.0
    control_states[:ev_volt_int] = 0.0
    
    # MPPT states
    control_states[:mppt_V_prev] = components.pv_system.Vmp
    control_states[:mppt_P_prev] = components.pv_system.Vmp * components.pv_system.Imp
    control_states[:mppt_integrator] = 0.0
    
    # VSG states
    control_states[:vsg_omega_int] = 0.0
    control_states[:vsg_voltage_int] = 0.0
    
    # AVR states
    control_states[:avr_integrator] = 0.0
    
    # PSS states
    control_states[:pss_states] = [0.0, 0.0, 0.0]
    
    return control_states
end

"""
Perform load flow calculation for initial conditions

Parameters:
- components: system components
- load_data: load specifications

Returns:
- bus_voltages: initial bus voltage phasors
- line_flows: initial line power flows
"""
function initial_load_flow(components, load_data=nothing)
    println("    - Performing initial load flow calculation...")
    
    # Simple 3-bus system for demonstration
    # Bus 1: Generator bus (slack)
    # Bus 2: Load bus
    # Bus 3: IBR bus
    
    # Bus data [bus_number, type, P_load, Q_load, V_mag, V_angle]
    # Type: 1=PQ, 2=PV, 3=Slack
    bus_data = [
        1  3  0.0   0.0   1.00  0.0;    # Slack bus (generator)
        2  1  0.6   0.3   1.00  0.0;    # Load bus
        3  1  0.2   0.1   1.00  0.0     # IBR bus
    ]
    
    # Line data [from_bus, to_bus, R, X, B/2]
    line_data = [
        1  2  0.01  0.10  0.05;
        2  3  0.02  0.15  0.03;
        1  3  0.015 0.12  0.04
    ]
    
    n_bus = size(bus_data, 1)
    
    # Build admittance matrix
    Y_bus = zeros(ComplexF64, n_bus, n_bus)
    
    for i in 1:size(line_data, 1)
        from_bus = Int(line_data[i, 1])
        to_bus = Int(line_data[i, 2])
        R = line_data[i, 3]
        X = line_data[i, 4]
        B = line_data[i, 5]
        
        Z = R + im * X
        y = 1.0 / Z
        
        # Add series admittance
        Y_bus[from_bus, to_bus] -= y
        Y_bus[to_bus, from_bus] -= y
        Y_bus[from_bus, from_bus] += y
        Y_bus[to_bus, to_bus] += y
        
        # Add shunt admittance
        Y_bus[from_bus, from_bus] += im * B
        Y_bus[to_bus, to_bus] += im * B
    end
    
    # Initial voltage guess
    V = ones(ComplexF64, n_bus)
    
    # Newton-Raphson load flow (simplified)
    max_iter = 20
    tolerance = 1e-6
    
    for iter in 1:max_iter
        # Calculate power mismatches
        S_calc = V .* conj.(Y_bus * V)
        P_calc = real.(S_calc)
        Q_calc = imag.(S_calc)
        
        # Power mismatches for PQ buses
        ΔP = zeros(n_bus-1)  # Exclude slack bus
        ΔQ = zeros(n_bus-1)
        
        bus_idx = 1
        for i in 2:n_bus  # Skip slack bus
            P_spec = -bus_data[i, 3]  # Load is negative injection
            Q_spec = -bus_data[i, 4]
            
            ΔP[bus_idx] = P_spec - P_calc[i]
            ΔQ[bus_idx] = Q_spec - Q_calc[i]
            bus_idx += 1
        end
        
        # Check convergence
        if maximum(abs.([ΔP; ΔQ])) < tolerance
            println("    - Load flow converged in $iter iterations")
            break
        end
        
        # Simplified voltage update (not full Newton-Raphson)
        for i in 2:n_bus
            V[i] *= (1.0 + 0.01 * (ΔP[i-1] + im * ΔQ[i-1]))
        end
    end
    
    # Calculate line flows
    line_flows = Dict()
    for i in 1:size(line_data, 1)
        from_bus = Int(line_data[i, 1])
        to_bus = Int(line_data[i, 2])
        
        R = line_data[i, 3]
        X = line_data[i, 4]
        Z = R + im * X
        
        I_flow = (V[from_bus] - V[to_bus]) / Z
        S_flow = V[from_bus] * conj(I_flow)
        
        line_flows["$from_bus-$to_bus"] = S_flow
    end
    
    bus_voltages = V
    
    println("    - Load flow calculation completed")
    println("    - Bus voltages: ", abs.(V))
    
    return bus_voltages, line_flows
end

"""
Initialize measurement and monitoring systems

Returns:
- measurement_config: configuration for measurement systems
"""
function initialize_measurements()
    measurement_config = Dict(
        # Voltage measurements
        :voltage_buses => [1, 2, 3],  # Buses with voltage measurements
        :voltage_noise => 0.001,      # Measurement noise (pu)
        
        # Current measurements  
        :current_lines => ["1-2", "2-3", "1-3"],  # Lines with current measurements
        :current_noise => 0.002,      # Measurement noise (pu)
        
        # Power measurements
        :power_generators => [1],     # Generators with power measurements
        :power_loads => [2, 3],       # Loads with power measurements
        :power_noise => 0.005,        # Measurement noise (pu)
        
        # Frequency measurements
        :frequency_buses => [1],      # Buses with frequency measurements
        :frequency_noise => 0.001,    # Measurement noise (Hz)
        
        # PMU measurements
        :pmu_buses => [1, 2],         # Buses with PMUs
        :pmu_sample_rate => 60,       # Samples per second
        :pmu_noise => 0.0005,         # PMU measurement noise
        
        # SCADA measurements
        :scada_sample_rate => 1,      # Samples per second
        :scada_delay => 2.0,          # Communication delay (s)
        
        # Protection measurements
        :protection_zones => ["zone1", "zone2", "zone3"],
        :relay_settings => Dict(
            "overcurrent" => 1.5,     # Pickup setting (pu)
            "undervoltage" => 0.8,    # Pickup setting (pu)
            "overfrequency" => 61.0,  # Pickup setting (Hz)
            "underfrequency" => 59.0  # Pickup setting (Hz)
        )
    )
    
    return measurement_config
end

"""
Initialize data logging and recording systems

Returns:
- logging_config: configuration for data logging
"""
function initialize_logging()
    logging_config = Dict(
        # Time series logging
        :log_states => true,
        :log_interval => DT * 10,     # Log every 10 time steps
        :log_variables => [
            :bus_voltages,
            :line_currents,
            :generator_power,
            :pv_power,
            :bess_power,
            :bess_soc,
            :motor_speed,
            :ev_soc,
            :system_frequency
        ],
        
        # Event logging
        :log_events => true,
        :event_types => [
            :fault_inception,
            :fault_clearing,
            :protection_operation,
            :control_action,
            :limit_violation
        ],
        
        # Disturbance recording
        :disturbance_recording => true,
        :trigger_conditions => [
            :voltage_deviation => 0.1,    # pu
            :frequency_deviation => 0.5,  # Hz
            :current_increase => 2.0      # pu
        ],
        :recording_duration => 10.0,     # seconds
        :pre_trigger_time => 1.0,        # seconds
        
        # File output settings
        :output_format => :hdf5,
        :output_directory => "./simulation_results/",
        :file_prefix => "em_transient_",
        :compression => true
    )
    
    return logging_config
end

"""
Validate initial conditions for physical consistency

Parameters:
- x0: initial state vector
- components: system components

Returns:
- validation_result: boolean indicating if validation passed
- validation_messages: array of validation messages
"""
function validate_initial_conditions(x0::Vector{Float64}, components)
    validation_result = true
    validation_messages = String[]
    
    println("  - Validating initial conditions...")
    
    # Check generator states
    gen_id = x0[STATE_INDICES[:gen_id]]
    gen_iq = x0[STATE_INDICES[:gen_iq]]
    gen_omega = x0[STATE_INDICES[:gen_ωr]]
    gen_delta = x0[STATE_INDICES[:gen_δ]]
    
    if abs(gen_omega - 1.0) > 0.1
        validation_result = false
        push!(validation_messages, "Generator speed deviation too large: $(gen_omega)")
    end
    
    if abs(gen_delta) > π/2
        validation_result = false
        push!(validation_messages, "Generator power angle too large: $(gen_delta)")
    end
    
    # Check PV states
    pv_Vdc = x0[STATE_INDICES[:pv_Vdc]]
    if pv_Vdc < 0 || pv_Vdc > 1000
        validation_result = false
        push!(validation_messages, "PV DC voltage out of range: $(pv_Vdc)")
    end
    
    # Check BESS states
    bess_SOC = x0[STATE_INDICES[:bess_SOC]]
    if bess_SOC < 0 || bess_SOC > 1
        validation_result = false
        push!(validation_messages, "BESS SOC out of range: $(bess_SOC)")
    end
    
    bess_Vdc = x0[STATE_INDICES[:bess_Vdc]]
    if bess_Vdc < 0 || bess_Vdc > 1200
        validation_result = false
        push!(validation_messages, "BESS DC voltage out of range: $(bess_Vdc)")
    end
    
    # Check motor states
    motor_omega = x0[STATE_INDICES[:motor_ωr]]
    if motor_omega < 0 || motor_omega > 2 * ω_BASE
        validation_result = false
        push!(validation_messages, "Motor speed out of range: $(motor_omega)")
    end
    
    # Check EV states
    ev_SOC = x0[STATE_INDICES[:ev_SOC]]
    if ev_SOC < 0 || ev_SOC > 1
        validation_result = false
        push!(validation_messages, "EV SOC out of range: $(ev_SOC)")
    end
    
    # Check network states
    bus_va = x0[STATE_INDICES[:bus_va]]
    bus_vb = x0[STATE_INDICES[:bus_vb]]
    bus_vc = x0[STATE_INDICES[:bus_vc]]
    
    V_mag = sqrt(bus_va^2 + bus_vb^2 + bus_vc^2) / √3
    if V_mag < 0.5 || V_mag > 1.5
        validation_result = false
        push!(validation_messages, "Bus voltage magnitude out of range: $(V_mag)")
    end
    
    # Check for NaN or Inf values
    for i in 1:length(x0)
        if !isfinite(x0[i])
            validation_result = false
            push!(validation_messages, "State x0[$i] is not finite: $(x0[i])")
        end
    end
    
    # Energy conservation check
    E_kinetic = 0.5 * components.generator.H * (gen_omega - 1.0)^2
    E_magnetic = 0.5 * components.line.L * sum(x0[STATE_INDICES[:line_ia]:STATE_INDICES[:line_ic]].^2)
    E_electric = 0.5 * components.line.C * sum(x0[STATE_INDICES[:bus_va]:STATE_INDICES[:bus_vc]].^2)
    E_battery = bess_SOC * components.bess.Q_max * components.bess.E0 / 1000  # kWh to MJ conversion
    
    E_total = E_kinetic + E_magnetic + E_electric + E_battery
    
    if E_total < 0
        validation_result = false
        push!(validation_messages, "Total system energy is negative: $(E_total)")
    end
    
    # Power balance check
    P_gen = gen_id * bus_va + gen_iq * bus_vb  # Simplified
    P_motor = x0[STATE_INDICES[:motor_ids]] * x0[STATE_INDICES[:motor_iqs]] * 0.1
    P_ev = x0[STATE_INDICES[:ev_iac]] * 240 / 1000  # kW
    P_balance = abs(P_gen - P_motor - P_ev - 0.1)  # 0.1 for other loads
    
    if P_balance > 0.2  # 20% tolerance
        validation_result = false
        push!(validation_messages, "Power balance error too large: $(P_balance)")
    end
    
    if validation_result
        println("  - Initial conditions validation passed!")
    else
        println("  - Initial conditions validation failed!")
        for msg in validation_messages
            println("    WARNING: $msg")
        end
    end
    
    return validation_result, validation_messages
end