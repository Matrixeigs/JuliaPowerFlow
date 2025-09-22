# components.jl - Component Structure Definitions

# ============================================================================
# COMPONENT STRUCTURES
# ============================================================================

# AC Generator (Synchronous Machine) Structure
struct SyncGenerator
    # Electrical parameters (in pu)
    Rs::Float64         # Stator resistance
    Xd::Float64         # d-axis synchronous reactance
    Xq::Float64         # q-axis synchronous reactance
    Xd_p::Float64       # d-axis transient reactance
    Xq_p::Float64       # q-axis transient reactance
    Xd_pp::Float64      # d-axis subtransient reactance
    Xq_pp::Float64      # q-axis subtransient reactance
    Xl::Float64         # Leakage reactance
    
    # Time constants
    Td0_p::Float64      # d-axis transient open-circuit time constant
    Tq0_p::Float64      # q-axis transient open-circuit time constant
    Td0_pp::Float64     # d-axis subtransient time constant
    Tq0_pp::Float64     # q-axis subtransient time constant
    
    # Mechanical parameters
    H::Float64          # Inertia constant (s)
    D::Float64          # Damping coefficient
    
    # Operating point
    V_t::Float64        # Terminal voltage (pu)
    P_gen::Float64      # Generated power (pu)
    Q_gen::Float64      # Generated reactive power (pu)
    
    # Constructor with default values
    function SyncGenerator(;Rs=GEN_PARAMS[:Rs], Xd=GEN_PARAMS[:Xd], 
                          Xq=GEN_PARAMS[:Xq], Xd_p=GEN_PARAMS[:Xd_p],
                          Xq_p=GEN_PARAMS[:Xq_p], Xd_pp=GEN_PARAMS[:Xd_pp],
                          Xq_pp=GEN_PARAMS[:Xq_pp], Xl=GEN_PARAMS[:Xl],
                          Td0_p=GEN_PARAMS[:Td0_p], Tq0_p=GEN_PARAMS[:Tq0_p],
                          Td0_pp=GEN_PARAMS[:Td0_pp], Tq0_pp=GEN_PARAMS[:Tq0_pp],
                          H=GEN_PARAMS[:H], D=GEN_PARAMS[:D],
                          V_t=GEN_PARAMS[:V_t], P_gen=GEN_PARAMS[:P_gen],
                          Q_gen=GEN_PARAMS[:Q_gen])
        new(Rs, Xd, Xq, Xd_p, Xq_p, Xd_pp, Xq_pp, Xl,
            Td0_p, Tq0_p, Td0_pp, Tq0_pp, H, D, V_t, P_gen, Q_gen)
    end
end

# PV System Structure
struct PVSystem
    # PV array parameters
    Isc::Float64        # Short circuit current (A)
    Voc::Float64        # Open circuit voltage (V)
    Imp::Float64        # Current at MPP (A)
    Vmp::Float64        # Voltage at MPP (V)
    Ki::Float64         # Temperature coefficient of current
    Kv::Float64         # Temperature coefficient of voltage
    Ns::Int             # Number of cells in series
    Np::Int             # Number of parallel strings
    
    # DC-DC converter parameters
    L_boost::Float64    # Boost inductor (H)
    C_dc::Float64       # DC capacitor (F)
    
    # Inverter parameters
    L_f::Float64        # Filter inductance (H)
    C_f::Float64        # Filter capacitance (F)
    R_f::Float64        # Filter resistance (Ω)
    
    # Control parameters
    Kp_mppt::Float64    # MPPT proportional gain
    Ki_mppt::Float64    # MPPT integral gain
    Kp_curr::Float64    # Current controller Kp
    Ki_curr::Float64    # Current controller Ki
    Kp_volt::Float64    # Voltage controller Kp
    Ki_volt::Float64    # Voltage controller Ki
    
    # Constructor with default values
    function PVSystem(;Isc=PV_PARAMS[:Isc], Voc=PV_PARAMS[:Voc],
                     Imp=PV_PARAMS[:Imp], Vmp=PV_PARAMS[:Vmp],
                     Ki=PV_PARAMS[:Ki], Kv=PV_PARAMS[:Kv],
                     Ns=PV_PARAMS[:Ns], Np=PV_PARAMS[:Np],
                     L_boost=PV_PARAMS[:L_boost], C_dc=PV_PARAMS[:C_dc],
                     L_f=PV_PARAMS[:L_f], C_f=PV_PARAMS[:C_f], R_f=PV_PARAMS[:R_f],
                     Kp_mppt=PV_PARAMS[:Kp_mppt], Ki_mppt=PV_PARAMS[:Ki_mppt],
                     Kp_curr=PV_PARAMS[:Kp_curr], Ki_curr=PV_PARAMS[:Ki_curr],
                     Kp_volt=PV_PARAMS[:Kp_volt], Ki_volt=PV_PARAMS[:Ki_volt])
        new(Isc, Voc, Imp, Vmp, Ki, Kv, Ns, Np, L_boost, C_dc,
            L_f, C_f, R_f, Kp_mppt, Ki_mppt, Kp_curr, Ki_curr, Kp_volt, Ki_volt)
    end
end

# BESS Structure
struct BESSSystem
    # Battery parameters
    E0::Float64         # No-load voltage (V)
    R_bat::Float64      # Internal resistance (Ω)
    K::Float64          # Polarization constant
    A::Float64          # Exponential zone amplitude
    B::Float64          # Exponential zone time constant
    Q_max::Float64      # Maximum capacity (Ah)
    
    # DC-DC converter parameters
    L_conv::Float64     # Converter inductance (H)
    C_dc::Float64       # DC link capacitance (F)
    
    # Inverter parameters
    L_f::Float64        # Filter inductance (H)
    C_f::Float64        # Filter capacitance (F)
    R_f::Float64        # Filter resistance (Ω)
    
    # Control parameters
    Kp_curr::Float64    # Current controller Kp
    Ki_curr::Float64    # Current controller Ki
    Kp_volt::Float64    # Voltage controller Kp
    Ki_volt::Float64    # Voltage controller Ki
    
    # Constructor with default values
    function BESSSystem(;E0=BESS_PARAMS[:E0], R_bat=BESS_PARAMS[:R_bat],
                       K=BESS_PARAMS[:K], A=BESS_PARAMS[:A], B=BESS_PARAMS[:B],
                       Q_max=BESS_PARAMS[:Q_max], L_conv=BESS_PARAMS[:L_conv],
                       C_dc=BESS_PARAMS[:C_dc], L_f=BESS_PARAMS[:L_f],
                       C_f=BESS_PARAMS[:C_f], R_f=BESS_PARAMS[:R_f],
                       Kp_curr=BESS_PARAMS[:Kp_curr], Ki_curr=BESS_PARAMS[:Ki_curr],
                       Kp_volt=BESS_PARAMS[:Kp_volt], Ki_volt=BESS_PARAMS[:Ki_volt])
        new(E0, R_bat, K, A, B, Q_max, L_conv, C_dc, L_f, C_f, R_f,
            Kp_curr, Ki_curr, Kp_volt, Ki_volt)
    end
end

# Transmission Line Structure
struct TransmissionLine
    length::Float64     # Line length (km)
    R::Float64         # Resistance per unit length (Ω/km)
    L::Float64         # Inductance per unit length (H/km)
    C::Float64         # Capacitance per unit length (F/km)
    G::Float64         # Conductance per unit length (S/km)
    
    # Constructor with default values
    function TransmissionLine(;length=LINE_PARAMS[:length], R=LINE_PARAMS[:R],
                             L=LINE_PARAMS[:L], C=LINE_PARAMS[:C], G=LINE_PARAMS[:G])
        new(length, R, L, C, G)
    end
end

# Transformer Structure
struct Transformer
    S_rated::Float64    # Rated power (VA)
    V1_rated::Float64   # Primary voltage (V)
    V2_rated::Float64   # Secondary voltage (V)
    R1::Float64        # Primary resistance (pu)
    R2::Float64        # Secondary resistance (pu)
    X1::Float64        # Primary leakage reactance (pu)
    X2::Float64        # Secondary leakage reactance (pu)
    Xm::Float64        # Magnetizing reactance (pu)
    Rm::Float64        # Core loss resistance (pu)
    
    # Constructor with default values
    function Transformer(;S_rated=TRAFO_PARAMS[:S_rated], V1_rated=TRAFO_PARAMS[:V1_rated],
                        V2_rated=TRAFO_PARAMS[:V2_rated], R1=TRAFO_PARAMS[:R1],
                        R2=TRAFO_PARAMS[:R2], X1=TRAFO_PARAMS[:X1],
                        X2=TRAFO_PARAMS[:X2], Xm=TRAFO_PARAMS[:Xm], Rm=TRAFO_PARAMS[:Rm])
        new(S_rated, V1_rated, V2_rated, R1, R2, X1, X2, Xm, Rm)
    end
end

# Induction Motor Structure
struct InductionMotor
    P_rated::Float64    # Rated power (W)
    V_rated::Float64    # Rated voltage (V)
    Rs::Float64        # Stator resistance (pu)
    Rr::Float64        # Rotor resistance (pu)
    Xs::Float64        # Stator leakage reactance (pu)
    Xr::Float64        # Rotor leakage reactance (pu)
    Xm::Float64        # Magnetizing reactance (pu)
    J::Float64         # Inertia (kg⋅m²)
    B::Float64         # Friction coefficient
    P::Int             # Number of poles
    
    # Constructor with default values
    function InductionMotor(;P_rated=MOTOR_PARAMS[:P_rated], V_rated=MOTOR_PARAMS[:V_rated],
                           Rs=MOTOR_PARAMS[:Rs], Rr=MOTOR_PARAMS[:Rr],
                           Xs=MOTOR_PARAMS[:Xs], Xr=MOTOR_PARAMS[:Xr],
                           Xm=MOTOR_PARAMS[:Xm], J=MOTOR_PARAMS[:J],
                           B=MOTOR_PARAMS[:B], P=MOTOR_PARAMS[:P])
        new(P_rated, V_rated, Rs, Rr, Xs, Xr, Xm, J, B, P)
    end
end

# EV Charger Structure
struct EVCharger
    # Battery parameters
    E0::Float64         # No-load voltage (V)
    R_bat::Float64      # Internal resistance (Ω)
    Q_max::Float64      # Maximum capacity (Ah)
    
    # Charger parameters
    P_max::Float64      # Maximum charging power (W)
    L_ac::Float64       # AC side inductance (H)
    C_dc::Float64       # DC link capacitance (F)
    
    # Control parameters
    Kp_curr::Float64    # Current controller Kp
    Ki_curr::Float64    # Current controller Ki
    
    # Constructor with default values
    function EVCharger(;E0=EV_PARAMS[:E0], R_bat=EV_PARAMS[:R_bat],
                      Q_max=EV_PARAMS[:Q_max], P_max=EV_PARAMS[:P_max],
                      L_ac=EV_PARAMS[:L_ac], C_dc=EV_PARAMS[:C_dc],
                      Kp_curr=EV_PARAMS[:Kp_curr], Ki_curr=EV_PARAMS[:Ki_curr])
        new(E0, R_bat, Q_max, P_max, L_ac, C_dc, Kp_curr, Ki_curr)
    end
end

# System Components Container
struct SystemComponents
    generator::SyncGenerator
    pv_system::PVSystem
    bess::BESSSystem
    line::TransmissionLine
    transformer::Transformer
    motor::InductionMotor
    ev_charger::EVCharger
    
    # Constructor
    function SystemComponents(gen, pv, bess, line, trafo, motor, ev)
        new(gen, pv, bess, line, trafo, motor, ev)
    end
end

# Environmental Conditions Structure
mutable struct EnvironmentalConditions
    irradiance::Float64     # Solar irradiance (W/m²)
    temperature::Float64    # Temperature (K)
    wind_speed::Float64     # Wind speed (m/s)
    
    function EnvironmentalConditions(;G=800.0, T=298.15, wind=5.0)
        new(G, T, wind)
    end
end

# Fault Condition Structure
mutable struct FaultCondition
    is_active::Bool
    fault_type::Symbol      # :three_phase, :line_to_ground, :line_to_line
    resistance::Float64     # Fault resistance (Ω)
    location::String        # Fault location description
    
    function FaultCondition(;active=false, type=:three_phase, R=0.001, loc="Bus1")
        new(active, type, R, loc)
    end
end