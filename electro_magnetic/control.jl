# control.jl - Control Systems for Inverters and Components

# ============================================================================
# CONTROL SYSTEM IMPLEMENTATIONS
# ============================================================================

"""
Grid-Forming Control using Virtual Synchronous Generator (VSG)

Parameters:
- P_meas: measured active power (pu)
- Q_meas: measured reactive power (pu)
- P_ref: reference active power (pu)
- Q_ref: reference reactive power (pu)
- ω_ref: reference frequency (pu)
- V_ref: reference voltage (pu)
- ω_int: frequency integrator state
- V_int: voltage integrator state
- dt: time step (s)

Returns:
- ω_new: new frequency reference (pu)
- V_new: new voltage reference (pu)
- ω_int_new: updated frequency integrator
- V_int_new: updated voltage integrator
"""
function grid_forming_vsg_control(P_meas::Float64, Q_meas::Float64, P_ref::Float64, Q_ref::Float64,
                                 ω_ref::Float64, V_ref::Float64, ω_int::Float64, V_int::Float64, dt::Float64)
    # VSG parameters
    J_vsg = 0.1      # Virtual inertia (s)
    D_vsg = 10.0     # Virtual damping
    m_p = CONTROL_PARAMS[:m_p]  # Active power droop
    m_q = CONTROL_PARAMS[:m_q]  # Reactive power droop
    
    # Active power control with virtual inertia and damping
    P_error = P_ref - P_meas
    ω_droop = ω_ref - m_p * P_error
    
    # Virtual swing equation
    dω_dt = (1/J_vsg) * (P_error - D_vsg * (ω_droop - ω_ref))
    ω_new = ω_droop + ω_int * dt
    ω_int_new = ω_int + dω_dt * dt
    
    # Reactive power control with voltage droop
    Q_error = Q_ref - Q_meas
    V_new = V_ref - m_q * Q_error + V_int
    V_int_new = V_int + 0.1 * Q_error * dt  # Small integral gain
    
    # Limit outputs
    ω_new = max(0.95, min(ω_new, 1.05))  # ±5% frequency deviation
    V_new = max(0.9, min(V_new, 1.1))    # ±10% voltage deviation
    
    return ω_new, V_new, ω_int_new, V_int_new
end

"""
Current Controller (PI) in DQ frame

Parameters:
- i_d_ref, i_q_ref: reference currents (A)
- i_d_meas, i_q_meas: measured currents (A)
- int_d, int_q: integrator states
- Kp, Ki: controller gains
- dt: time step (s)
- ω: angular frequency (rad/s)
- L: inductance (H)

Returns:
- v_d_ref, v_q_ref: reference voltages (V)
- int_d_new, int_q_new: updated integrator states
"""
function current_controller_dq(i_d_ref::Float64, i_q_ref::Float64, i_d_meas::Float64, i_q_meas::Float64,
                              int_d::Float64, int_q::Float64, Kp::Float64, Ki::Float64, dt::Float64,
                              ω::Float64, L::Float64)
    # Current errors
    e_d = i_d_ref - i_d_meas
    e_q = i_q_ref - i_q_meas
    
    # PI controller with decoupling
    v_d_ref = Kp * e_d + Ki * int_d - ω * L * i_q_meas  # Cross-coupling compensation
    v_q_ref = Kp * e_q + Ki * int_q + ω * L * i_d_meas  # Cross-coupling compensation
    
    # Update integrators with anti-windup
    int_d_new = int_d + e_d * dt
    int_q_new = int_q + e_q * dt
    
    # Anti-windup (simple clamping)
    v_max = 1000.0  # Maximum voltage
    if abs(v_d_ref) > v_max
        int_d_new = int_d  # Stop integrator
        v_d_ref = sign(v_d_ref) * v_max
    end
    if abs(v_q_ref) > v_max
        int_q_new = int_q  # Stop integrator
        v_q_ref = sign(v_q_ref) * v_max
    end
    
    return v_d_ref, v_q_ref, int_d_new, int_q_new
end

"""
Voltage Controller (PI) for Grid-Forming Inverters

Parameters:
- v_d_ref, v_q_ref: reference voltages (V)
- v_d_meas, v_q_meas: measured voltages (V)
- int_d, int_q: integrator states
- Kp, Ki: controller gains
- dt: time step (s)

Returns:
- i_d_ref, i_q_ref: reference currents (A)
- int_d_new, int_q_new: updated integrator states
"""
function voltage_controller_dq(v_d_ref::Float64, v_q_ref::Float64, v_d_meas::Float64, v_q_meas::Float64,
                              int_d::Float64, int_q::Float64, Kp::Float64, Ki::Float64, dt::Float64)
    # Voltage errors
    e_d = v_d_ref - v_d_meas
    e_q = v_q_ref - v_q_meas
    
    # PI controller
    i_d_ref = Kp * e_d + Ki * int_d
    i_q_ref = Kp * e_q + Ki * int_q
    
    # Update integrators
    int_d_new = int_d + e_d * dt
    int_q_new = int_q + e_q * dt
    
    # Current limiting
    i_max = 100.0  # Maximum current (A)
    if abs(i_d_ref) > i_max
        int_d_new = int_d  # Stop integrator
        i_d_ref = sign(i_d_ref) * i_max
    end
    if abs(i_q_ref) > i_max
        int_q_new = int_q  # Stop integrator
        i_q_ref = sign(i_q_ref) * i_max
    end
    
    return i_d_ref, i_q_ref, int_d_new, int_q_new
end

"""
Phase-Locked Loop (PLL) for Grid Synchronization

Parameters:
- v_abc: three-phase voltage measurements [Va, Vb, Vc]
- θ_pll: PLL angle state (rad)
- ω_pll: PLL frequency state (rad/s)
- int_pll: PLL integrator state
- Kp, Ki: PLL gains
- dt: time step (s)

Returns:
- θ_pll_new: updated PLL angle (rad)
- ω_pll_new: updated PLL frequency (rad/s)
- int_pll_new: updated integrator state
"""
function phase_locked_loop(v_abc::Vector{Float64}, θ_pll::Float64, ω_pll::Float64, int_pll::Float64,
                          Kp::Float64, Ki::Float64, dt::Float64)
    # Transform to dq frame
    v_dq0 = abc_to_dq0(v_abc, θ_pll)
    v_d, v_q = v_dq0[1], v_dq0[2]
    
    # PLL error (q-component should be zero when locked)
    error = v_q
    
    # PI controller
    ω_pll_new = ω_BASE + Kp * error + Ki * int_pll
    
    # Update integrator
    int_pll_new = int_pll + error * dt
    
    # Update angle
    θ_pll_new = θ_pll + ω_pll_new * dt
    
    # Wrap angle to [0, 2π]
    θ_pll_new = mod(θ_pll_new, 2π)
    
    # Frequency limiting
    ω_pll_new = max(0.9 * ω_BASE, min(ω_pll_new, 1.1 * ω_BASE))
    
    return θ_pll_new, ω_pll_new, int_pll_new
end

"""
DC Voltage Controller for DC-DC Converters

Parameters:
- V_dc_ref: reference DC voltage (V)
- V_dc_meas: measured DC voltage (V)
- int_dc: integrator state
- Kp, Ki: controller gains
- dt: time step (s)

Returns:
- duty_cycle: PWM duty cycle (0-1)
- int_dc_new: updated integrator state
"""
function dc_voltage_controller(V_dc_ref::Float64, V_dc_meas::Float64, int_dc::Float64,
                              Kp::Float64, Ki::Float64, dt::Float64)
    # Voltage error
    error = V_dc_ref - V_dc_meas
    
    # PI controller
    duty_cycle = Kp * error + Ki * int_dc
    
    # Update integrator
    int_dc_new = int_dc + error * dt
    
    # Duty cycle limiting with anti-windup
    if duty_cycle > 0.95
        duty_cycle = 0.95
        int_dc_new = int_dc  # Stop integrator
    elseif duty_cycle < 0.05
        duty_cycle = 0.05
        int_dc_new = int_dc  # Stop integrator
    end
    
    return duty_cycle, int_dc_new
end

"""
Battery State of Charge (SOC) Controller

Parameters:
- SOC_ref: reference SOC (0-1)
- SOC_meas: measured SOC (0-1)
- I_max: maximum charging/discharging current (A)
- int_soc: integrator state
- Kp, Ki: controller gains
- dt: time step (s)

Returns:
- I_ref: reference current (A) - positive for discharge
- int_soc_new: updated integrator state
"""
function soc_controller(SOC_ref::Float64, SOC_meas::Float64, I_max::Float64, int_soc::Float64,
                       Kp::Float64, Ki::Float64, dt::Float64)
    # SOC error
    error = SOC_ref - SOC_meas
    
    # PI controller
    I_ref = Kp * error + Ki * int_soc
    
    # Update integrator
    int_soc_new = int_soc + error * dt
    
    # Current limiting with anti-windup
    if I_ref > I_max
        I_ref = I_max
        int_soc_new = int_soc  # Stop integrator
    elseif I_ref < -I_max
        I_ref = -I_max
        int_soc_new = int_soc  # Stop integrator
    end
    
    # SOC-based current limiting (protect battery)
    if SOC_meas > 0.95 && I_ref < 0  # Nearly full, limit charging
        I_ref = max(I_ref, -0.1 * I_max)
    elseif SOC_meas < 0.05 && I_ref > 0  # Nearly empty, limit discharging
        I_ref = min(I_ref, 0.1 * I_max)
    end
    
    return I_ref, int_soc_new
end

"""
Motor Speed Controller

Parameters:
- ω_ref: reference speed (rad/s)
- ω_meas: measured speed (rad/s)
- int_speed: integrator state
- Kp, Ki: controller gains
- dt: time step (s)

Returns:
- T_ref: reference torque (Nm)
- int_speed_new: updated integrator state
"""
function motor_speed_controller(ω_ref::Float64, ω_meas::Float64, int_speed::Float64,
                               Kp::Float64, Ki::Float64, dt::Float64)
    # Speed error
    error = ω_ref - ω_meas
    
    # PI controller
    T_ref = Kp * error + Ki * int_speed
    
    # Update integrator
    int_speed_new = int_speed + error * dt
    
    # Torque limiting
    T_max = 1000.0  # Maximum torque (Nm)
    if abs(T_ref) > T_max
        T_ref = sign(T_ref) * T_max
        int_speed_new = int_speed  # Stop integrator
    end
    
    return T_ref, int_speed_new
end

"""
Automatic Voltage Regulator (AVR) for Synchronous Generator

Parameters:
- V_ref: reference terminal voltage (pu)
- V_meas: measured terminal voltage (pu)
- int_avr: integrator state
- Kp, Ki: controller gains
- dt: time step (s)

Returns:
- V_fd: field voltage (pu)
- int_avr_new: updated integrator state
"""
function automatic_voltage_regulator(V_ref::Float64, V_meas::Float64, int_avr::Float64,
                                    Kp::Float64, Ki::Float64, dt::Float64)
    # Voltage error
    error = V_ref - V_meas
    
    # PI controller
    V_fd = Kp * error + Ki * int_avr
    
    # Update integrator
    int_avr_new = int_avr + error * dt
    
    # Field voltage limiting
    V_fd_max = 5.0  # Maximum field voltage (pu)
    V_fd_min = 0.0  # Minimum field voltage (pu)
    
    if V_fd > V_fd_max
        V_fd = V_fd_max
        int_avr_new = int_avr  # Stop integrator
    elseif V_fd < V_fd_min
        V_fd = V_fd_min
        int_avr_new = int_avr  # Stop integrator
    end
    
    return V_fd, int_avr_new
end

"""
Power System Stabilizer (PSS) for Damping Enhancement

Parameters:
- ω_meas: measured rotor speed (pu)
- ω_ref: reference speed (pu)
- pss_states: PSS internal states [x1, x2, x3]
- dt: time step (s)

Returns:
- V_pss: PSS output voltage (pu)
- pss_states_new: updated PSS states
"""
function power_system_stabilizer(ω_meas::Float64, ω_ref::Float64, pss_states::Vector{Float64}, dt::Float64)
    # PSS parameters
    K_pss = 20.0    # PSS gain
    T1 = 0.05       # Lead time constant
    T2 = 0.02       # Lag time constant
    T3 = 3.0        # Washout time constant
    
    # Speed deviation
    Δω = ω_meas - ω_ref
    
    # Washout filter (high-pass)
    x1 = pss_states[1]
    dx1_dt = (Δω - x1) / T3
    x1_new = x1 + dx1_dt * dt
    
    # Lead-lag compensator
    x2, x3 = pss_states[2], pss_states[3]
    
    # First lead-lag stage
    dx2_dt = (x1_new - x2) / T2
    x2_new = x2 + dx2_dt * dt
    y1 = x2_new + T1 * dx2_dt
    
    # Second lead-lag stage (if needed)
    dx3_dt = (y1 - x3) / T2
    x3_new = x3 + dx3_dt * dt
    V_pss = K_pss * (x3_new + T1 * dx3_dt)
    
    # Output limiting
    V_pss_max = 0.1  # ±10% of rated voltage
    V_pss = max(-V_pss_max, min(V_pss, V_pss_max))
    
    pss_states_new = [x1_new, x2_new, x3_new]
    
    return V_pss, pss_states_new
end

"""
Fault Detection and Protection Logic

Parameters:
- I_abc: three-phase current measurements [Ia, Ib, Ic]
- V_abc: three-phase voltage measurements [Va, Vb, Vc]
- I_threshold: overcurrent threshold (A)
- V_threshold: undervoltage threshold (V)

Returns:
- fault_detected: boolean indicating fault detection
- fault_type: type of fault detected
- protection_action: recommended protection action
"""
function fault_detection_protection(I_abc::Vector{Float64}, V_abc::Vector{Float64},
                                   I_threshold::Float64, V_threshold::Float64)
    # Calculate RMS values
    I_rms = [abs(I) for I in I_abc]
    V_rms = [abs(V) for V in V_abc]
    
    # Overcurrent detection
    overcurrent = any(I_rms .> I_threshold)
    
    # Undervoltage detection
    undervoltage = any(V_rms .< V_threshold)
    
    # Unbalance detection
    I_avg = mean(I_rms)
    V_avg = mean(V_rms)
    I_unbalance = maximum(abs.(I_rms .- I_avg)) / I_avg > 0.1
    V_unbalance = maximum(abs.(V_rms .- V_avg)) / V_avg > 0.02
    
    # Determine fault type and action
    fault_detected = false
    fault_type = :none
    protection_action = :none
    
    if overcurrent
        fault_detected = true
        if I_rms[1] > I_threshold && I_rms[2] > I_threshold && I_rms[3] > I_threshold
            fault_type = :three_phase_fault
            protection_action = :trip_immediately
        elseif sum(I_rms .> I_threshold) == 2
            fault_type = :line_to_line_fault
            protection_action = :trip_delayed
        else
            fault_type = :single_phase_fault
            protection_action = :trip_delayed
        end
    elseif undervoltage
        fault_detected = true
        fault_type = :undervoltage
        protection_action = :disconnect_load
    elseif I_unbalance || V_unbalance
        fault_detected = true
        fault_type = :unbalance
        protection_action = :alarm_only
    end
    
    return fault_detected, fault_type, protection_action
end