# parameters.jl - System Parameters and Constants

# ============================================================================
# GLOBAL SYSTEM PARAMETERS
# ============================================================================

# Base values
const V_BASE = 13800.0      # Base voltage (V)
const S_BASE = 10e6         # Base power (VA)
const F_BASE = 60.0         # Base frequency (Hz)
const ω_BASE = 2π * F_BASE  # Base angular frequency
const Z_BASE = V_BASE^2 / S_BASE
const I_BASE = S_BASE / (√3 * V_BASE)

# Time parameters
const DT = 1e-5      # Time step (10 μs for electromagnetic transients)
const T_SIM = 0.5    # Simulation time (s)

# Environmental parameters
const G_REF = 1000.0    # Reference irradiance (W/m²)
const T_REF = 298.15    # Reference temperature (K)

# Physical constants
const Q_ELECTRON = 1.602176634e-19  # Elementary charge (C)
const K_BOLTZMANN = 1.380649e-23    # Boltzmann constant (J/K)

# ============================================================================
# COMPONENT-SPECIFIC PARAMETERS
# ============================================================================

# Synchronous Generator Default Parameters
const GEN_PARAMS = Dict(
    :Rs => 0.003,      # Stator resistance (pu)
    :Xd => 1.8,        # d-axis synchronous reactance (pu)
    :Xq => 1.7,        # q-axis synchronous reactance (pu)
    :Xd_p => 0.3,      # d-axis transient reactance (pu)
    :Xq_p => 0.55,     # q-axis transient reactance (pu)
    :Xd_pp => 0.25,    # d-axis subtransient reactance (pu)
    :Xq_pp => 0.25,    # q-axis subtransient reactance (pu)
    :Xl => 0.2,        # Leakage reactance (pu)
    :Td0_p => 8.0,     # d-axis transient time constant (s)
    :Tq0_p => 0.4,     # q-axis transient time constant (s)
    :Td0_pp => 0.03,   # d-axis subtransient time constant (s)
    :Tq0_pp => 0.05,   # q-axis subtransient time constant (s)
    :H => 3.7,         # Inertia constant (s)
    :D => 0.0,         # Damping coefficient
    :V_t => 1.0,       # Terminal voltage (pu)
    :P_gen => 0.8,     # Generated power (pu)
    :Q_gen => 0.3      # Generated reactive power (pu)
)

# PV System Default Parameters
const PV_PARAMS = Dict(
    :Isc => 8.21,          # Short circuit current (A)
    :Voc => 32.9,          # Open circuit voltage (V)
    :Imp => 7.61,          # Current at MPP (A)
    :Vmp => 26.3,          # Voltage at MPP (V)
    :Ki => 0.0032,         # Temperature coefficient of current
    :Kv => -0.123,         # Temperature coefficient of voltage
    :Ns => 54,             # Number of cells in series
    :Np => 1,              # Number of parallel strings
    :L_boost => 1e-3,      # Boost inductor (H)
    :C_dc => 2200e-6,      # DC capacitor (F)
    :L_f => 1e-3,          # Filter inductance (H)
    :C_f => 10e-6,         # Filter capacitance (F)
    :R_f => 0.1,           # Filter resistance (Ω)
    :Kp_mppt => 0.1,       # MPPT proportional gain
    :Ki_mppt => 10.0,      # MPPT integral gain
    :Kp_curr => 10.0,      # Current controller Kp
    :Ki_curr => 1000.0,    # Current controller Ki
    :Kp_volt => 1.0,       # Voltage controller Kp
    :Ki_volt => 100.0      # Voltage controller Ki
)

# BESS Default Parameters
const BESS_PARAMS = Dict(
    :E0 => 12.6,           # No-load voltage (V)
    :R_bat => 0.01,        # Internal resistance (Ω)
    :K => 0.1,             # Polarization constant
    :A => 0.3,             # Exponential zone amplitude
    :B => 3.0,             # Exponential zone time constant
    :Q_max => 100.0,       # Maximum capacity (Ah)
    :L_conv => 1e-3,       # Converter inductance (H)
    :C_dc => 4400e-6,      # DC link capacitance (F)
    :L_f => 1e-3,          # Filter inductance (H)
    :C_f => 10e-6,         # Filter capacitance (F)
    :R_f => 0.1,           # Filter resistance (Ω)
    :Kp_curr => 10.0,      # Current controller Kp
    :Ki_curr => 1000.0,    # Current controller Ki
    :Kp_volt => 1.0,       # Voltage controller Kp
    :Ki_volt => 100.0      # Voltage controller Ki
)

# Transmission Line Default Parameters
const LINE_PARAMS = Dict(
    :length => 10.0,       # Line length (km)
    :R => 0.1,             # Resistance per unit length (Ω/km)
    :L => 1e-3,            # Inductance per unit length (H/km)
    :C => 10e-9,           # Capacitance per unit length (F/km)
    :G => 1e-6             # Conductance per unit length (S/km)
)

# Transformer Default Parameters
const TRAFO_PARAMS = Dict(
    :S_rated => 10e6,      # Rated power (VA)
    :V1_rated => 13800,    # Primary voltage (V)
    :V2_rated => 480,      # Secondary voltage (V)
    :R1 => 0.01,           # Primary resistance (pu)
    :R2 => 0.01,           # Secondary resistance (pu)
    :X1 => 0.08,           # Primary leakage reactance (pu)
    :X2 => 0.08,           # Secondary leakage reactance (pu)
    :Xm => 50.0,           # Magnetizing reactance (pu)
    :Rm => 500.0           # Core loss resistance (pu)
)

# Induction Motor Default Parameters
const MOTOR_PARAMS = Dict(
    :P_rated => 100e3,     # Rated power (W)
    :V_rated => 480,       # Rated voltage (V)
    :Rs => 0.01,           # Stator resistance (pu)
    :Rr => 0.007,          # Rotor resistance (pu)
    :Xs => 0.05,           # Stator leakage reactance (pu)
    :Xr => 0.05,           # Rotor leakage reactance (pu)
    :Xm => 3.0,            # Magnetizing reactance (pu)
    :J => 2.0,             # Inertia (kg⋅m²)
    :B => 0.01,            # Friction coefficient
    :P => 4                # Number of poles
)

# EV Charger Default Parameters
const EV_PARAMS = Dict(
    :E0 => 400.0,          # No-load voltage (V)
    :R_bat => 0.1,         # Internal resistance (Ω)
    :Q_max => 75.0,        # Maximum capacity (Ah)
    :P_max => 50e3,        # Maximum charging power (W)
    :L_ac => 1e-3,         # AC side inductance (H)
    :C_dc => 1000e-6,      # DC link capacitance (F)
    :Kp_curr => 5.0,       # Current controller Kp
    :Ki_curr => 500.0      # Current controller Ki
)

# Control System Parameters
const CONTROL_PARAMS = Dict(
    :m_p => 0.05,          # Active power droop coefficient
    :m_q => 0.05,          # Reactive power droop coefficient
    :τ => 1.0,             # Integration time constant
    :ω_ref => 1.0,         # Reference frequency (pu)
    :V_ref => 1.0          # Reference voltage (pu)
)

# Fault Parameters
const FAULT_PARAMS = Dict(
    :t_start => 0.1,       # Fault start time (s)
    :t_end => 0.15,        # Fault end time (s)
    :R_fault => 0.001,     # Fault resistance (Ω)
    :fault_type => :three_phase  # Fault type
)

# State Vector Indices (for easy reference)
const STATE_INDICES = Dict(
    # Generator states (7 states)
    :gen_id => 1,
    :gen_iq => 2,
    :gen_ifd => 3,
    :gen_ikd => 4,
    :gen_ikq => 5,
    :gen_ωr => 6,
    :gen_δ => 7,
    
    # PV states (5 states)
    :pv_iL => 8,
    :pv_Vdc => 9,
    :pv_mppt_int => 10,
    :pv_curr_int_d => 11,
    :pv_curr_int_q => 12,
    
    # BESS states (8 states)
    :bess_ibat => 13,
    :bess_SOC => 14,
    :bess_Vdc => 15,
    :bess_curr_int_d => 16,
    :bess_curr_int_q => 17,
    :bess_volt_int_d => 18,
    :bess_volt_int_q => 19,
    :bess_ω_int => 20,
    
    # Motor states (5 states)
    :motor_ids => 21,
    :motor_iqs => 22,
    :motor_idr => 23,
    :motor_iqr => 24,
    :motor_ωr => 25,
    
    # EV states (5 states)
    :ev_iac => 26,
    :ev_SOC => 27,
    :ev_Vdc => 28,
    :ev_curr_int => 29,
    :ev_volt_int => 30,
    
    # Network states (10 states)
    :line_ia => 31,
    :line_ib => 32,
    :line_ic => 33,
    :bus_va => 34,
    :bus_vb => 35,
    :bus_vc => 36,
    :trafo_flux_a => 37,
    :trafo_flux_b => 38,
    :trafo_flux_c => 39,
    :fault_current => 40
)

const N_STATES = 40  # Total number of states