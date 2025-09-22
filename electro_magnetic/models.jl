# models.jl - Mathematical Models for System Components

# ============================================================================
# COMPONENT MATHEMATICAL MODELS
# ============================================================================

"""
PV Cell I-V Characteristic using Single Diode Model

Parameters:
- V_pv: PV voltage (V)
- G: irradiance (W/m²)
- T: temperature (K)
- pv: PV system parameters

Returns:
- I_pv: PV current (A)
"""
function pv_current_model(V_pv::Float64, G::Float64, T::Float64, pv::PVSystem)
    # Physical constants
    q = Q_ELECTRON
    k = K_BOLTZMANN
    
    # Reference conditions
    G_ref = G_REF
    T_ref = T_REF
    
    # Temperature-dependent parameters
    I_ph = pv.Np * (pv.Isc + pv.Ki * (T - T_ref)) * G / G_ref
    V_t = k * T / q
    a = 1.3  # Ideality factor
    
    # Reverse saturation current
    I_0 = pv.Np * pv.Isc / (exp(pv.Voc / (a * V_t * pv.Ns)) - 1)
    
    # Series and shunt resistance (estimated from datasheet)
    R_s = (pv.Voc - pv.Vmp) / pv.Imp
    R_sh = pv.Vmp / (pv.Isc - pv.Imp)
    
    # Newton-Raphson solution for I_pv
    function current_equation(I)
        return I_ph - I_0 * (exp((V_pv + I * R_s / pv.Np) / (a * V_t * pv.Ns)) - 1) - 
               (V_pv + I * R_s / pv.Np) / (R_sh * pv.Ns / pv.Np) - I
    end
    
    function current_derivative(I)
        exp_term = exp((V_pv + I * R_s / pv.Np) / (a * V_t * pv.Ns))
        return -I_0 * (R_s / pv.Np) / (a * V_t * pv.Ns) * exp_term - 
               (R_s / pv.Np) / (R_sh * pv.Ns / pv.Np) - 1
    end
    
    # Initial guess
    I_pv = I_ph * 0.9
    
    # Newton-Raphson iterations
    for i in 1:20
        f_val = current_equation(I_pv)
        if abs(f_val) < 1e-8
            break
        end
        df_val = current_derivative(I_pv)
        if abs(df_val) < 1e-12
            break
        end
        I_new = I_pv - f_val / df_val
        I_pv = max(0.0, min(I_new, I_ph))  # Constrain to physical limits
    end
    
    return max(0.0, I_pv)
end

"""
Battery Voltage Model using Shepherd Model

Parameters:
- I_bat: battery current (A) - positive for discharge
- SOC: state of charge (0-1)
- Q: extracted capacity (Ah)
- bess: BESS parameters

Returns:
- V_bat: battery voltage (V)
"""
function battery_voltage_model(I_bat::Float64, SOC::Float64, Q::Float64, bess::BESSSystem)
    if SOC <= 0.0
        return 0.0
    end
    
    # Polarization resistance (simplified)
    K_polarization = bess.K * bess.Q_max / max(bess.Q_max - Q, 0.1)
    
    # Exponential zone voltage
    V_exp = bess.A * exp(-bess.B * Q)
    
    # Battery voltage equation
    V_bat = bess.E0 - bess.R_bat * I_bat - K_polarization * I_bat + V_exp
    
    # Apply SOC-dependent voltage correction
    V_soc_correction = 0.1 * (SOC - 0.5)  # Simple linear correction
    V_bat += V_soc_correction
    
    return max(0.0, V_bat)
end

"""
Synchronous Generator Electrical Model (Park's Equations)

Parameters:
- x: state vector [id, iq, ifd, ikd, ikq, ωr, δ]
- gen: generator parameters

Returns:
- derivatives: [did/dt, diq/dt, difd/dt, dikd/dt, dikq/dt]
"""
function sync_generator_electrical_model(x::Vector{Float64}, gen::SyncGenerator)
    id, iq, ifd, ikd, ikq, ωr, δ = x[1:7]
    
    # Flux linkages
    ψd = -gen.Xd * id + gen.Xd * ifd - gen.Xd * ikd
    ψq = -gen.Xq * iq - gen.Xq * ikq
    ψfd = gen.Xd * ifd + gen.Xd * id - gen.Xd * ikd
    ψkd = -gen.Xd * ikd + gen.Xd * ifd + gen.Xd * id
    ψkq = -gen.Xq * ikq - gen.Xq * iq
    
    # Terminal voltage (assumed constant for grid connection)
    vd = 0.0  # Assuming infinite bus
    vq = 1.0  # Per unit terminal voltage
    
    # Stator equations
    did_dt = (1/gen.Xd) * (-gen.Rs * id - ωr * ψq + vd)
    diq_dt = (1/gen.Xq) * (-gen.Rs * iq + ωr * ψd + vq)
    
    # Field winding equation
    vfd = 1.0  # Field voltage (constant excitation)
    difd_dt = (1/gen.Td0_p) * (vfd - gen.Rs * ifd)
    
    # Damper winding equations
    dikd_dt = -(1/gen.Td0_pp) * ikd
    dikq_dt = -(1/gen.Tq0_pp) * ikq
    
    return [did_dt, diq_dt, difd_dt, dikd_dt, dikq_dt]
end

"""
Synchronous Generator Mechanical Model

Parameters:
- x: state vector [id, iq, ifd, ikd, ikq, ωr, δ]
- gen: generator parameters
- T_m: mechanical torque (pu)

Returns:
- [dωr/dt, dδ/dt]
"""
function sync_generator_mechanical_model(x::Vector{Float64}, gen::SyncGenerator, T_m::Float64)
    id, iq, ωr, δ = x[1], x[2], x[6], x[7]
    
    # Electromagnetic torque
    T_e = id * iq  # Simplified expression
    
    # Swing equation
    dωr_dt = (1/(2*gen.H)) * (T_m - T_e - gen.D * (ωr - 1.0))
    dδ_dt = ω_BASE * (ωr - 1.0)
    
    return [dωr_dt, dδ_dt]
end

"""
Induction Motor Model (dq frame)

Parameters:
- x: state vector [ids, iqs, idr, iqr, ωr]
- motor: motor parameters
- vds, vqs: stator voltages
- T_load: load torque

Returns:
- derivatives: [dids/dt, diqs/dt, didr/dt, diqr/dt, dωr/dt]
"""
function induction_motor_model(x::Vector{Float64}, motor::InductionMotor, 
                              vds::Float64, vqs::Float64, T_load::Float64)
    ids, iqs, idr, iqr, ωr = x
    
    # Flux linkages
    ψds = motor.Xs * ids + motor.Xm * idr
    ψqs = motor.Xs * iqs + motor.Xm * iqr
    ψdr = motor.Xr * idr + motor.Xm * ids
    ψqr = motor.Xr * iqr + motor.Xm * iqs
    
    # Slip frequency
    ωs = ω_BASE  # Synchronous frequency
    ωsl = ωs - ωr  # Slip frequency
    
    # Stator equations
    dids_dt = (1/motor.Xs) * (vds - motor.Rs * ids + ωs * ψqs)
    diqs_dt = (1/motor.Xs) * (vqs - motor.Rs * iqs - ωs * ψds)
    
    # Rotor equations (short-circuited)
    didr_dt = (1/motor.Xr) * (-motor.Rr * idr + ωsl * ψqr)
    diqr_dt = (1/motor.Xr) * (-motor.Rr * iqr - ωsl * ψdr)
    
    # Electromagnetic torque
    T_e = (3/2) * (motor.P/2) * (ψds * iqs - ψqs * ids)
    
    # Mechanical equation
    dωr_dt = (1/motor.J) * (T_e - T_load - motor.B * ωr)
    
    return [dids_dt, diqs_dt, didr_dt, diqr_dt, dωr_dt]
end

"""
Transmission Line Model (Bergeron's Method)

Parameters:
- I_send: sending end current
- I_recv: receiving end current
- V_send: sending end voltage
- V_recv: receiving end voltage
- line: line parameters

Returns:
- line equations for distributed parameter model
"""
function transmission_line_model(I_send::Float64, I_recv::Float64, 
                                V_send::Float64, V_recv::Float64, line::TransmissionLine)
    # Line parameters
    R = line.R * line.length
    L = line.L * line.length
    C = line.C * line.length
    G = line.G * line.length
    
    # Characteristic impedance (resistive approximation)
    Z_c = sqrt((R + im*ω_BASE*L)/(G + im*ω_BASE*C))
    Z_c_real = real(Z_c)
    
    # Propagation constant
    γ = sqrt((R + im*ω_BASE*L)*(G + im*ω_BASE*C))
    
    # Travel time
    τ = line.length / 3e8  # Speed of light approximation
    
    # Bergeron equivalent circuit (simplified)
    # This is a simplified representation - full implementation would require
    # convolution integrals and history terms
    
    dI_send_dt = (1/L) * (V_send - R * I_send - V_recv)
    dI_recv_dt = (1/L) * (V_recv - R * I_recv - V_send)
    
    return [dI_send_dt, dI_recv_dt]
end

"""
Transformer Model with Saturation

Parameters:
- i_p: primary current vector [ia, ib, ic]
- i_s: secondary current vector [ia, ib, ic]
- flux: magnetizing flux vector [φa, φb, φc]
- trafo: transformer parameters

Returns:
- transformer differential equations
"""
function transformer_model(i_p::Vector{Float64}, i_s::Vector{Float64}, 
                         flux::Vector{Float64}, trafo::Transformer)
    # Turns ratio
    n = trafo.V1_rated / trafo.V2_rated
    
    # Magnetizing current (with saturation)
    i_mag = saturation_curve.(flux, trafo.Xm)
    
    # Primary equations
    v_p = zeros(3)  # Would be provided by network solution
    di_p_dt = (1/trafo.X1) .* (v_p - trafo.R1 .* i_p - flux)
    
    # Secondary equations  
    v_s = zeros(3)  # Would be provided by network solution
    di_s_dt = (1/trafo.X2) .* (v_s - trafo.R2 .* i_s - flux./n)
    
    # Flux equations
    dflux_dt = v_p - trafo.R1 .* i_p - trafo.X1 .* di_p_dt
    
    return [di_p_dt; di_s_dt; dflux_dt]
end

"""
Saturation Curve for Transformer Core

Parameters:
- flux: magnetic flux (Wb)
- X_m: magnetizing reactance (unsaturated)

Returns:
- i_mag: magnetizing current (A)
"""
function saturation_curve(flux::Float64, X_m::Float64)
    # Simple exponential saturation model
    flux_knee = 1.0  # Knee point flux (pu)
    sat_factor = 2.0  # Saturation factor
    
    if abs(flux) <= flux_knee
        i_mag = flux / X_m
    else
        i_mag_linear = flux_knee / X_m
        flux_excess = abs(flux) - flux_knee
        i_mag_sat = i_mag_linear + flux_excess / (X_m / sat_factor)
        i_mag = sign(flux) * i_mag_sat
    end
    
    return i_mag
end

"""
Environmental Conditions Update

Parameters:
- t: current time (s)

Returns:
- env: updated environmental conditions
"""
function update_environmental_conditions(t::Float64)
    # Time-varying irradiance (cloud effects)
    G_base = 800.0
    G_variation = 200.0 * sin(0.1 * t) + 100.0 * sin(0.5 * t)
    G = max(0.0, G_base + G_variation)
    
    # Time-varying temperature
    T_base = 298.15  # 25°C
    T_variation = 10.0 * sin(0.05 * t)
    T = T_base + T_variation
    
    # Wind speed (affects cooling)
    wind_base = 5.0
    wind_variation = 3.0 * sin(0.2 * t)
    wind = max(0.0, wind_base + wind_variation)
    
    return EnvironmentalConditions(G=G, T=T, wind=wind)
end

"""
MPPT Control Algorithm (Perturb and Observe)

Parameters:
- V_pv: PV voltage (V)
- I_pv: PV current (A)
- V_pv_prev: previous PV voltage (V)
- P_pv_prev: previous PV power (W)
- step_size: voltage step size (V)

Returns:
- V_ref_new: new voltage reference (V)
- P_pv: current PV power (W)
"""
function mppt_perturb_observe(V_pv::Float64, I_pv::Float64, V_pv_prev::Float64, 
                             P_pv_prev::Float64, step_size::Float64)
    P_pv = V_pv * I_pv
    
    if P_pv > P_pv_prev
        if V_pv > V_pv_prev
            V_ref_new = V_pv + step_size  # Move in same direction
        else
            V_ref_new = V_pv - step_size  # Move in same direction
        end
    else
        if V_pv > V_pv_prev
            V_ref_new = V_pv - step_size  # Reverse direction
        else
            V_ref_new = V_pv + step_size  # Reverse direction
        end
    end
    
    # Limit voltage reference
    V_ref_new = max(0.1, min(V_ref_new, 50.0))
    
    return V_ref_new, P_pv
end