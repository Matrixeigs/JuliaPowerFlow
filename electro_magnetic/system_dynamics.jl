# system_dynamics.jl - Complete System Differential Equations

# ============================================================================
# MAIN SYSTEM DYNAMICS FUNCTION
# ============================================================================

"""
Complete System Differential Equations

This function implements the full electromagnetic transient model for the
AC distribution system with inverter-based resources.

Parameters:
- dx: derivative vector (output)
- x: state vector (input)
- p: parameter tuple containing all system components
- t: current time (s)

State Vector Organization:
x[1:7]   - Generator states [id, iq, ifd, ikd, ikq, ωr, δ]
x[8:12]  - PV states [iL, Vdc, mppt_int, curr_int_d, curr_int_q]
x[13:20] - BESS states [ibat, SOC, Vdc, curr_int_d, curr_int_q, volt_int_d, volt_int_q, ω_int]
x[21:25] - Motor states [ids, iqs, idr, iqr, ωr]
x[26:30] - EV states [iac, SOC, Vdc, curr_int, volt_int]
x[31:40] - Network states [line currents, bus voltages, transformer fluxes, fault current]
"""
function system_dynamics!(dx, x, p, t)
    # Unpack system components
    components = p
    gen = components.generator
    pv = components.pv_system
    bess = components.bess
    line = components.line
    trafo = components.transformer
    motor = components.motor
    ev = components.ev_charger
    
    # Initialize derivative vector
    fill!(dx, 0.0)
    
    # Update environmental conditions
    env = update_environmental_conditions(t)
    
    # Check for fault conditions
    fault = FaultCondition()
    if FAULT_PARAMS[:t_start] <= t <= FAULT_PARAMS[:t_end]
        fault.is_active = true
        fault.fault_type = FAULT_PARAMS[:fault_type]
        fault.resistance = FAULT_PARAMS[:R_fault]
    end
    
    # ========================================================================
    # SYNCHRONOUS GENERATOR DYNAMICS
    # ========================================================================
    
    # Extract generator states
    gen_states = x[STATE_INDICES[:gen_id]:STATE_INDICES[:gen_δ]]
    
    # Electrical dynamics
    gen_elec_derivs = sync_generator_electrical_model(gen_states, gen)
    dx[STATE_INDICES[:gen_id]:STATE_INDICES[:gen_ikq]] = gen_elec_derivs
    
    # Mechanical dynamics
    T_m = 0.8  # Mechanical torque (constant prime mover)
    gen_mech_derivs = sync_generator_mechanical_model(gen_states, gen, T_m)
    dx[STATE_INDICES[:gen_ωr]:STATE_INDICES[:gen_δ]] = gen_mech_derivs
    
    # ========================================================================
    # PV SYSTEM DYNAMICS
    # ========================================================================
    
    # Extract PV states
    iL_pv = x[STATE_INDICES[:pv_iL]]
    Vdc_pv = x[STATE_INDICES[:pv_Vdc]]
    mppt_int = x[STATE_INDICES[:pv_mppt_int]]
    curr_int_d_pv = x[STATE_INDICES[:pv_curr_int_d]]
    curr_int_q_pv = x[STATE_INDICES[:pv_curr_int_q]]
    
    # PV array model
    V_pv = 26.0  # Simplified PV voltage (would be calculated from MPPT)
    I_pv = pv_current_model(V_pv, env.irradiance, env.temperature, pv)
    
    # DC-DC boost converter dynamics
    d_boost = 0.6  # Duty cycle from MPPT controller
    dx[STATE_INDICES[:pv_iL]] = (1/pv.L_boost) * (V_pv - (1 - d_boost) * Vdc_pv)
    
    # DC link dynamics
    I_inv_pv = 5.0  # Inverter input current (simplified)
    dx[STATE_INDICES[:pv_Vdc]] = (1/pv.C_dc) * ((1 - d_boost) * iL_pv - I_inv_pv)
    
    # MPPT integrator
    P_pv = V_pv * I_pv
    P_ref_pv = 200.0  # Reference power
    dx[STATE_INDICES[:pv_mppt_int]] = P_ref_pv - P_pv
    
    # Current controller integrators (simplified)
    dx[STATE_INDICES[:pv_curr_int_d]] = 0.0  # Would be error from current controller
    dx[STATE_INDICES[:pv_curr_int_q]] = 0.0
    
    # ========================================================================
    # BESS DYNAMICS
    # ========================================================================
    
    # Extract BESS states
    ibat = x[STATE_INDICES[:bess_ibat]]
    SOC_bess = x[STATE_INDICES[:bess_SOC]]
    Vdc_bess = x[STATE_INDICES[:bess_Vdc]]
    curr_int_d_bess = x[STATE_INDICES[:bess_curr_int_d]]
    curr_int_q_bess = x[STATE_INDICES[:bess_curr_int_q]]
    volt_int_d_bess = x[STATE_INDICES[:bess_volt_int_d]]
    volt_int_q_bess = x[STATE_INDICES[:bess_volt_int_q]]
    ω_int_bess = x[STATE_INDICES[:bess_ω_int]]
    
    # Battery model
    Q_extracted = bess.Q_max * (1 - SOC_bess)
    V_bat = battery_voltage_model(ibat, SOC_bess, Q_extracted, bess)
    
    # DC-DC converter dynamics (bidirectional)
    dx[STATE_INDICES[:bess_ibat]] = (1/bess.L_conv) * (V_bat - Vdc_bess)
    
    # SOC dynamics
    dx[STATE_INDICES[:bess_SOC]] = -ibat / bess.Q_max
    
    # DC link dynamics
    I_inv_bess = 3.0  # Inverter current (simplified)
    dx[STATE_INDICES[:bess_Vdc]] = (1/bess.C_dc) * (ibat - I_inv_bess)
    
    # Grid-forming control integrators
    P_meas_bess = V_bat * ibat
    Q_meas_bess = 0.0  # Simplified
    P_ref_bess = 100.0
    Q_ref_bess = 0.0
    
    ω_new, V_new, ω_int_new, V_int_new = grid_forming_vsg_control(
        P_meas_bess, Q_meas_bess, P_ref_bess, Q_ref_bess,
        1.0, 1.0, ω_int_bess, 0.0, DT
    )
    
    dx[STATE_INDICES[:bess_ω_int]] = ω_int_new - ω_int_bess
    
    # Current and voltage controller integrators (simplified)
    dx[STATE_INDICES[:bess_curr_int_d]] = 0.0
    dx[STATE_INDICES[:bess_curr_int_q]] = 0.0
    dx[STATE_INDICES[:bess_volt_int_d]] = 0.0
    dx[STATE_INDICES[:bess_volt_int_q]] = 0.0
    
    # ========================================================================
    # INDUCTION MOTOR DYNAMICS
    # ========================================================================
    
    # Extract motor states
    motor_states = x[STATE_INDICES[:motor_ids]:STATE_INDICES[:motor_ωr]]
    
    # Motor supply voltages (from network solution)
    vds_motor = 1.0  # Simplified - would come from network
    vqs_motor = 0.0
    T_load = 50.0    # Load torque (Nm)
    
    # Motor dynamics
    motor_derivs = induction_motor_model(motor_states, motor, vds_motor, vqs_motor, T_load)
    dx[STATE_INDICES[:motor_ids]:STATE_INDICES[:motor_ωr]] = motor_derivs
    
    # ========================================================================
    # EV CHARGER DYNAMICS
    # ========================================================================
    
    # Extract EV states
    iac_ev = x[STATE_INDICES[:ev_iac]]
    SOC_ev = x[STATE_INDICES[:ev_SOC]]
    Vdc_ev = x[STATE_INDICES[:ev_Vdc]]
    curr_int_ev = x[STATE_INDICES[:ev_curr_int]]
    volt_int_ev = x[STATE_INDICES[:ev_volt_int]]
    
    # EV battery model (simplified)
    V_bat_ev = ev.E0 - ev.R_bat * iac_ev
    
    # Charging dynamics
    if SOC_ev < 0.9  # Charge if not full
        I_ref_ev = 10.0  # Charging current reference
    else
        I_ref_ev = 0.0   # Stop charging when full
    end
    
    # AC side dynamics
    V_grid_ev = 240.0  # Grid voltage (V)
    dx[STATE_INDICES[:ev_iac]] = (1/ev.L_ac) * (V_grid_ev - 240.0 * iac_ev / I_BASE)
    
    # SOC dynamics
    dx[STATE_INDICES[:ev_SOC]] = iac_ev / ev.Q_max
    
    # DC link dynamics
    dx[STATE_INDICES[:ev_Vdc]] = (1/ev.C_dc) * (iac_ev - V_bat_ev / ev.E0)
    
    # system_dynamics.jl - Complete System Differential Equations (Continued)

    # Controller integrators
    dx[STATE_INDICES[:ev_curr_int]] = I_ref_ev - iac_ev
    dx[STATE_INDICES[:ev_volt_int]] = 400.0 - Vdc_ev
    
    # ========================================================================
    # NETWORK DYNAMICS
    # ========================================================================
    
    # Extract network states
    line_ia = x[STATE_INDICES[:line_ia]]
    line_ib = x[STATE_INDICES[:line_ib]]
    line_ic = x[STATE_INDICES[:line_ic]]
    bus_va = x[STATE_INDICES[:bus_va]]
    bus_vb = x[STATE_INDICES[:bus_vb]]
    bus_vc = x[STATE_INDICES[:bus_vc]]
    trafo_flux_a = x[STATE_INDICES[:trafo_flux_a]]
    trafo_flux_b = x[STATE_INDICES[:trafo_flux_b]]
    trafo_flux_c = x[STATE_INDICES[:trafo_flux_c]]
    fault_current = x[STATE_INDICES[:fault_current]]
    
    # Transmission line dynamics (simplified π-model)
    R_line = line.R * line.length
    L_line = line.L * line.length
    C_line = line.C * line.length
    
    # Line current dynamics (simplified)
    V_send = [bus_va, bus_vb, bus_vc]
    V_recv = [0.95*bus_va, 0.95*bus_vb, 0.95*bus_vc]  # Simplified receiving end
    I_line = [line_ia, line_ib, line_ic]
    
    for phase in 1:3
        phase_idx = STATE_INDICES[:line_ia] + phase - 1
        dx[phase_idx] = (1/L_line) * (V_send[phase] - V_recv[phase] - R_line * I_line[phase])
    end
    
    # Bus voltage dynamics (capacitive)
    I_inject_a = gen_states[1] + 0.1 * sin(ω_BASE * t)  # Generator + load injection
    I_inject_b = gen_states[1] + 0.1 * sin(ω_BASE * t - 2π/3)
    I_inject_c = gen_states[1] + 0.1 * sin(ω_BASE * t + 2π/3)
    
    dx[STATE_INDICES[:bus_va]] = (1/C_line) * (I_inject_a - line_ia)
    dx[STATE_INDICES[:bus_vb]] = (1/C_line) * (I_inject_b - line_ib)
    dx[STATE_INDICES[:bus_vc]] = (1/C_line) * (I_inject_c - line_ic)
    
    # Transformer flux dynamics
    V_primary = [bus_va, bus_vb, bus_vc]
    I_primary = [line_ia, line_ib, line_ic]
    
    for phase in 1:3
        flux_idx = STATE_INDICES[:trafo_flux_a] + phase - 1
        dx[flux_idx] = V_primary[phase] - trafo.R1 * I_primary[phase]
    end
    
    # ========================================================================
    # FAULT DYNAMICS
    # ========================================================================
    
    if fault.is_active
        # Three-phase short circuit dynamics
        if fault.fault_type == :three_phase
            # Fault current growth limited by system impedance
            Z_fault = fault.resistance + 0.1  # Total fault impedance
            V_fault = sqrt(bus_va^2 + bus_vb^2 + bus_vc^2) / √3
            
            dx[STATE_INDICES[:fault_current]] = (V_fault - Z_fault * fault_current) / 0.001
            
            # Fault affects bus voltages
            fault_effect = fault_current * 0.1
            dx[STATE_INDICES[:bus_va]] -= fault_effect
            dx[STATE_INDICES[:bus_vb]] -= fault_effect
            dx[STATE_INDICES[:bus_vc]] -= fault_effect
            
            # Fault affects generator (increased current)
            dx[STATE_INDICES[:gen_id]] += fault_current * 0.01
            dx[STATE_INDICES[:gen_iq]] += fault_current * 0.01
        end
    else
        # No fault - fault current decays
        dx[STATE_INDICES[:fault_current]] = -fault_current * 10.0
    end
    
    # ========================================================================
    # COUPLING AND INTERACTION EFFECTS
    # ========================================================================
    
    # Generator-network coupling
    # The generator terminal voltage affects the network
    V_gen_terminal = sqrt(gen_states[1]^2 + gen_states[2]^2)
    network_coupling = 0.1 * (V_gen_terminal - 1.0)
    dx[STATE_INDICES[:bus_va]] += network_coupling
    
    # PV-network coupling through inverter
    # PV inverter output affects network voltage
    P_pv_out = Vdc_pv * I_inv_pv * 0.95  # Inverter efficiency
    pv_voltage_effect = P_pv_out * 0.001
    dx[STATE_INDICES[:bus_va]] += pv_voltage_effect
    
    # BESS-network coupling (grid-forming)
    # BESS provides voltage support
    V_bess_ref = V_new  # From grid-forming control
    bess_voltage_support = 0.05 * (V_bess_ref - 1.0)
    dx[STATE_INDICES[:bus_va]] += bess_voltage_support
    dx[STATE_INDICES[:bus_vb]] += bess_voltage_support
    dx[STATE_INDICES[:bus_vc]] += bess_voltage_support
    
    # Motor load effect on network
    P_motor = motor_states[1] * motor_states[2] * motor.V_rated
    motor_load_effect = P_motor / (V_BASE * I_BASE) * 0.01
    dx[STATE_INDICES[:bus_va]] -= motor_load_effect
    dx[STATE_INDICES[:bus_vb]] -= motor_load_effect
    dx[STATE_INDICES[:bus_vc]] -= motor_load_effect
    
    # EV charging load effect
    P_ev_charging = Vdc_ev * iac_ev
    ev_load_effect = P_ev_charging / (V_BASE * I_BASE) * 0.005
    dx[STATE_INDICES[:bus_va]] -= ev_load_effect
    
    # ========================================================================
    # PROTECTION AND CONTROL ACTIONS
    # ========================================================================
    
    # Monitor system conditions for protection
    I_abc_monitor = [line_ia, line_ib, line_ic]
    V_abc_monitor = [bus_va, bus_vb, bus_vc]
    
    fault_detected, fault_type_detected, protection_action = fault_detection_protection(
        I_abc_monitor, V_abc_monitor, 2.0 * I_BASE, 0.8 * V_BASE
    )
    
    if fault_detected && protection_action == :trip_immediately
        # Emergency shutdown of inverters
        dx[STATE_INDICES[:pv_iL]] *= 0.1      # Rapid PV shutdown
        dx[STATE_INDICES[:bess_ibat]] *= 0.1   # BESS emergency stop
        dx[STATE_INDICES[:ev_iac]] *= 0.1      # EV charger disconnect
    end
    
    # ========================================================================
    # FREQUENCY AND VOLTAGE REGULATION
    # ========================================================================
    
    # System frequency calculation (simplified)
    ω_system = gen_states[6] * ω_BASE
    
    # Frequency-dependent load shedding
    if ω_system < 0.98 * ω_BASE
        # Under-frequency load shedding
        load_shed_factor = 0.9
        dx[STATE_INDICES[:motor_ids]] *= load_shed_factor
        dx[STATE_INDICES[:motor_iqs]] *= load_shed_factor
        dx[STATE_INDICES[:ev_iac]] *= load_shed_factor
    end
    
    # Voltage regulation through reactive power
    V_avg = (abs(bus_va) + abs(bus_vb) + abs(bus_vc)) / 3
    if V_avg < 0.95
        # Voltage support from BESS
        dx[STATE_INDICES[:bess_curr_int_q]] += 0.1 * (0.95 - V_avg)
    elseif V_avg > 1.05
        # Voltage reduction
        dx[STATE_INDICES[:bess_curr_int_q]] -= 0.1 * (V_avg - 1.05)
    end
    
    # ========================================================================
    # HARMONIC AND POWER QUALITY EFFECTS
    # ========================================================================
    
    # Add harmonic distortion from inverters (simplified)
    harmonic_freq = 5 * ω_BASE  # 5th harmonic
    harmonic_amplitude = 0.05   # 5% of fundamental
    
    # Harmonic injection from PV inverter
    pv_harmonic = harmonic_amplitude * sin(harmonic_freq * t)
    dx[STATE_INDICES[:bus_va]] += pv_harmonic * 0.01
    
    # Harmonic injection from BESS inverter
    bess_harmonic = harmonic_amplitude * sin(harmonic_freq * t + π/3)
    dx[STATE_INDICES[:bus_vb]] += bess_harmonic * 0.01
    
    # Interharmonic effects
    interharmonic_freq = 3.5 * ω_BASE
    interharmonic = 0.02 * sin(interharmonic_freq * t)
    dx[STATE_INDICES[:bus_vc]] += interharmonic * 0.005
    
    # ========================================================================
    # TEMPERATURE AND AGING EFFECTS
    # ========================================================================
    
    # Temperature effects on component parameters (simplified)
    temp_factor = 1.0 + 0.001 * (env.temperature - T_REF)
    
    # Battery aging effect on capacity
    if t > 100.0  # After 100 seconds of operation
        aging_factor = 0.999  # Slight capacity reduction
        dx[STATE_INDICES[:bess_SOC]] *= aging_factor
    end
    
    # PV panel temperature effects
    pv_temp_effect = 1.0 - 0.004 * (env.temperature - T_REF)
    dx[STATE_INDICES[:pv_iL]] *= pv_temp_effect
    
    # ========================================================================
    # STOCHASTIC DISTURBANCES
    # ========================================================================
    
    # Add small random disturbances to represent measurement noise and
    # unmodeled dynamics
    noise_level = 1e-4
    
    for i in 1:length(dx)
        if abs(dx[i]) > 1e-6  # Only add noise to active states
            dx[i] += noise_level * randn() * abs(dx[i])
        end
    end
    
    # ========================================================================
    # NUMERICAL STABILITY CHECKS
    # ========================================================================
    
    # Limit derivative magnitudes to prevent numerical instability
    max_derivative = 1e6
    for i in 1:length(dx)
        if abs(dx[i]) > max_derivative
            dx[i] = sign(dx[i]) * max_derivative
        end
    end
    
    # Check for NaN or Inf values
    for i in 1:length(dx)
        if !isfinite(dx[i])
            dx[i] = 0.0
        end
    end
    
    return nothing
end

"""
Algebraic Constraint Equations for Network Solution

This function implements the algebraic constraints that must be satisfied
at each time step, representing Kirchhoff's laws and power balance.

Parameters:
- x: current state vector
- t: current time
- components: system components

Returns:
- constraints: vector of constraint violations
"""
function network_constraints(x::Vector{Float64}, t::Float64, components)
    constraints = zeros(10)  # Number of algebraic constraints
    
    # Extract relevant states
    gen_id = x[STATE_INDICES[:gen_id]]
    gen_iq = x[STATE_INDICES[:gen_iq]]
    line_ia = x[STATE_INDICES[:line_ia]]
    line_ib = x[STATE_INDICES[:line_ib]]
    line_ic = x[STATE_INDICES[:line_ic]]
    bus_va = x[STATE_INDICES[:bus_va]]
    bus_vb = x[STATE_INDICES[:bus_vb]]
    bus_vc = x[STATE_INDICES[:bus_vc]]
    
    # Kirchhoff's Current Law at main bus
    I_gen_a = gen_id * cos(ω_BASE * t) - gen_iq * sin(ω_BASE * t)
    I_gen_b = gen_id * cos(ω_BASE * t - 2π/3) - gen_iq * sin(ω_BASE * t - 2π/3)
    I_gen_c = gen_id * cos(ω_BASE * t + 2π/3) - gen_iq * sin(ω_BASE * t + 2π/3)
    
    constraints[1] = I_gen_a - line_ia - 0.1  # Load current
    constraints[2] = I_gen_b - line_ib - 0.1
    constraints[3] = I_gen_c - line_ic - 0.1
    
    # Kirchhoff's Voltage Law
    constraints[4] = bus_va - 1.0 * cos(ω_BASE * t)      # Reference voltage
    constraints[5] = bus_vb - 1.0 * cos(ω_BASE * t - 2π/3)
    constraints[6] = bus_vc - 1.0 * cos(ω_BASE * t + 2π/3)
    
    # Power balance constraints
    P_gen = gen_id * bus_va + gen_iq * bus_vb  # Simplified
    P_load = 0.8  # Total system load
    constraints[7] = P_gen - P_load
    
    # Voltage magnitude constraints
    V_mag = sqrt(bus_va^2 + bus_vb^2 + bus_vc^2) / √3
    constraints[8] = V_mag - 1.0  # Per unit voltage
    
    # Frequency constraint
    ω_system = x[STATE_INDICES[:gen_ωr]] * ω_BASE
    constraints[9] = ω_system - ω_BASE
    
    # Phase balance constraint
    constraints[10] = bus_va + bus_vb + bus_vc  # Should sum to zero
    
    return constraints
end

"""
Event Detection Function for Discrete Events

This function detects events that require special handling during simulation,
such as fault inception, clearing, or protection operations.

Parameters:
- x: current state vector
- t: current time
- integrator: ODE integrator object

Returns:
- event_occurred: boolean indicating if event occurred
- event_type: type of event
"""
function detect_events(x::Vector{Float64}, t::Float64, integrator)
    event_occurred = false
    event_type = :none
    
    # Fault inception detection
    if abs(t - FAULT_PARAMS[:t_start]) < 1e-6
        event_occurred = true
        event_type = :fault_inception
    end
    
    # Fault clearing detection
    if abs(t - FAULT_PARAMS[:t_end]) < 1e-6
        event_occurred = true
        event_type = :fault_clearing
    end
    
    # Overcurrent detection
    I_total = sqrt(x[STATE_INDICES[:gen_id]]^2 + x[STATE_INDICES[:gen_iq]]^2)
    if I_total > 3.0  # 3 pu overcurrent
        event_occurred = true
        event_type = :overcurrent_trip
    end
    
    # Undervoltage detection
    V_bus = abs(x[STATE_INDICES[:bus_va]])
    if V_bus < 0.7  # 70% voltage
        event_occurred = true
        event_type = :undervoltage_trip
    end
    
    # BESS SOC limits
    SOC_bess = x[STATE_INDICES[:bess_SOC]]
    if SOC_bess < 0.05 || SOC_bess > 0.95
        event_occurred = true
        event_type = :bess_limit
    end
    
    return event_occurred, event_type
end