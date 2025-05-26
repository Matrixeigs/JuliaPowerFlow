# Induction Motor Mathematical Model

## Overview

This document presents the complete mathematical model for three-phase induction motors, including the equivalent circuit, electrical characteristics, mechanical relationships, and the characteristic U-P-Q relationships that are crucial for power system analysis.

## 1. Equivalent Circuit Model

### 1.1 Per-Phase Equivalent Circuit

The induction motor is represented by the following per-phase equivalent circuit:

```
V₁ ----R₁----jX₁----+----R₂'/s----jX₂'----
                     |
                    jXₘ
                     |
                    ---
```

Where:
- **V₁**: Stator phase voltage (V)
- **R₁**: Stator resistance per phase (Ω)
- **X₁**: Stator leakage reactance per phase (Ω)
- **R₂'**: Rotor resistance referred to stator per phase (Ω)
- **X₂'**: Rotor leakage reactance referred to stator per phase (Ω)
- **Xₘ**: Magnetizing reactance per phase (Ω)
- **s**: Slip

### 1.2 Complex Impedances

```julia
Z₁ = R₁ + jX₁                    # Stator impedance
Z₂ = R₂'/s + jX₂'               # Rotor impedance (referred to stator)
Zₘ = jXₘ                        # Magnetizing impedance
```

### 1.3 Total Circuit Impedance

The parallel combination of rotor and magnetizing branches:
```julia
Z_parallel = (Z₂ × Zₘ) / (Z₂ + Zₘ)
Z_total = Z₁ + Z_parallel
```

## 2. Electrical Characteristics

### 2.1 Current Calculations

**Stator Current:**
```julia
I₁ = V₁ / Z_total
```

**Voltage across parallel combination:**
```julia
V_parallel = I₁ × Z_parallel
```

**Rotor Current (referred to stator):**
```julia
I₂ = V_parallel / Z₂
```

**Magnetizing Current:**
```julia
Iₘ = V_parallel / Zₘ
```

**Current relationship verification:**
```julia
I₁ = I₂ + Iₘ    # Kirchhoff's current law
```

### 2.2 Power Analysis

**Three-Phase Power Calculations:**

**Stator Copper Loss:**
```julia
P₁_loss = 3 × |I₁|² × R₁
```

**Rotor Copper Loss:**
```julia
P₂_loss = 3 × |I₂|² × R₂'
```

**Air Gap Power:**
```julia
P_airgap = 3 × |I₂|² × R₂'/s
```

**Mechanical Power Output:**
```julia
P_mech = P_airgap - P₂_loss = P_airgap × (1 - s)
```

**Electrical Power Input:**
```julia
P_elec = P_airgap + P₁_loss
```

### 2.3 Reactive Power Analysis

The reactive power has three main components that create the characteristic U-shaped Q vs V curve:

**Stator Leakage Reactive Power:**
```julia
Q₁ = 3 × |I₁|² × X₁
```

**Rotor Leakage Reactive Power:**
```julia
Q₂ = 3 × |I₂|² × X₂'
```

**Magnetizing Reactive Power:**
```julia
Qₘ = 3 × |V_parallel|² / Xₘ
```

**Total Reactive Power:**
```julia
Q_total = Q₁ + Q₂ + Qₘ
```

**Alternative calculation using complex power:**
```julia
S_total = 3 × V₁ × I₁*
Q_total = Im(S_total)
```

## 3. Slip and Speed Relationships

### 3.1 Synchronous Speed

```julia
nₛ = 120 × f / p    # rpm
```
Where:
- **f**: Supply frequency (Hz)
- **p**: Number of poles

### 3.2 Slip Definition

```julia
s = (nₛ - n) / nₛ
```
Where:
- **n**: Rotor speed (rpm)
- **nₛ**: Synchronous speed (rpm)

### 3.3 Angular Velocities

```julia
ωₛ = 2π × nₛ / 60    # Synchronous angular velocity (rad/s)
ωᵣ = 2π × n / 60     # Rotor angular velocity (rad/s)
ωᵣ = ωₛ × (1 - s)    # Relationship between angular velocities
```

## 4. Mechanical Characteristics

### 4.1 Torque Calculation

**Electromagnetic Torque:**
```julia
T = P_airgap / ωₛ
```

**Alternative torque expression:**
```julia
T = (3 × |I₂|² × R₂'/s) / ωₛ
```

### 4.2 Torque-Slip Relationship

For analysis purposes, the torque-slip relationship can be approximated as:
```julia
T = (3 × V₁² × R₂'/s) / (ωₛ × [(R₁ + R₂'/s)² + (X₁ + X₂')²])
```

### 4.3 Maximum Torque

The maximum torque occurs at slip:
```julia
s_max_torque = R₂' / √(R₁² + (X₁ + X₂')²)
```

## 5. Terminal Voltage Effects (U-P-Q Relationship)

### 5.1 Voltage Dependency and Load Considerations

**Two Operating Scenarios:**

1. **Constant Slip (Theoretical Analysis):**
   - Slip remains fixed regardless of voltage
   - Useful for understanding circuit behavior
   - Not realistic for actual motor operation

2. **Constant Load (Realistic Operation):**
   - Motor maintains constant mechanical torque
   - Slip varies with voltage to maintain load
   - Represents real motor behavior

### 5.2 Constant Load Analysis (Realistic Behavior)

**Slip Variation with Voltage:**

For constant torque output, the required slip changes with voltage:

```julia
# At low voltage: Higher slip needed for same torque
s_required ∝ 1/V²    # Approximately

# This affects rotor resistance:
R₂_effective = R₂'/s_required

# At low voltage:
# s↑ → R₂'/s↓ → Higher rotor current → More losses
```

**Physical Explanation:**
- **Low Voltage**: To maintain torque, slip increases significantly
- **Higher Slip**: Reduces effective rotor resistance R₂'/s
- **Lower R₂'/s**: Increases rotor current for same power transfer
- **Higher Current**: Increases I²R losses in rotor
- **Result**: Efficiency drops dramatically at low voltages

**Efficiency with Constant Load:**
```julia
# Rotor losses increase with higher slip:
P_rotor_loss = 3 × I₂² × R₂'

# Where I₂ increases due to lower R₂'/s:
I₂ = V_parallel / (R₂'/s + jX₂')

# At low voltage: s↑ → R₂'/s↓ → I₂↑ → P_loss↑ → η↓
```

### 5.3 Comparison: Constant Slip vs Constant Load

**Constant Slip (Unrealistic):**
- Efficiency appears stable across voltage range
- R₂'/s remains constant
- Does not represent actual motor behavior

**Constant Load (Realistic):**
- Efficiency drops significantly at low voltage
- Slip increases dramatically below rated voltage
- Matches observed motor behavior in practice

**Key Relationships:**
```julia
# Constant slip scenario:
η_constant_slip ≈ constant (unrealistic)

# Constant load scenario:
η_constant_load = f(V²/s(V))  # Decreases with low voltage
```

### 5.4 Practical Voltage Effects

**Operating Recommendations:**
1. **Maintain rated voltage**: Critical for efficiency
2. **Avoid operation below 0.9 p.u.**: Efficiency drops rapidly
3. **Consider slip increase**: Low voltage requires higher slip
4. **Monitor rotor losses**: Increase significantly at low voltage

**Voltage Drop Effects:**
- **90% voltage**: ~20% efficiency reduction
- **80% voltage**: ~40% efficiency reduction  
- **70% voltage**: ~60% efficiency reduction (severe)

## 6. Efficiency Calculations

### 6.1 Motor Efficiency

```julia
η = P_mech / P_elec × 100%
```

### 6.2 Loss Components

**Total Losses:**
```julia
P_losses = P₁_loss + P₂_loss + P_core + P_friction + P_stray
```

Where:
- **P_core**: Core losses (approximately constant)
- **P_friction**: Friction and windage losses
- **P_stray**: Stray load losses

## 7. Implementation in Julia

### 7.1 Motor Structure Definition

```julia
mutable struct InductionMotor
    # Nameplate data
    rated_power::Float64      # kW
    rated_voltage::Float64    # V
    rated_frequency::Float64  # Hz
    rated_speed::Float64      # rpm
    poles::Int               # Number of poles
    
    # Equivalent circuit parameters
    r1::Float64              # Stator resistance (Ω)
    x1::Float64              # Stator reactance (Ω)
    r2::Float64              # Rotor resistance (Ω)
    x2::Float64              # Rotor reactance (Ω)
    xm::Float64              # Magnetizing reactance (Ω)
end
```

### 7.2 Performance Calculation Function with Saturation

```julia
function calculate_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)
    # Get saturated magnetizing reactance - CRITICAL for U-shape
    xm_saturated = calculate_saturated_xm(motor, voltage)
    
    # Equivalent circuit analysis with saturation
    V1 = voltage / √3
    Z1 = motor.r1 + im * motor.x1
    Z2 = motor.r2/slip + im * motor.x2
    Zm = im * xm_saturated  # Use saturated value
    
    # Circuit solution
    Z_parallel = (Z2 * Zm) / (Z2 + Zm)
    Z_total = Z1 + Z_parallel
    I1 = V1 / Z_total
    
    # Reactive power with saturation effects
    Q_stator = 3 * abs(I1)² * motor.x1
    Q_rotor = 3 * abs(I1 * Z_parallel / Z2)² * motor.x2  
    Q_magnetizing = 3 * abs(I1 * Z_parallel)² / xm_saturated
    
    Q_total = Q_stator + Q_rotor + Q_magnetizing
    
    return P_mech, Q_total, abs(I1), efficiency
end
```

### 7.3 Constant Load Analysis

```julia
function calculate_motor_slip_for_load(motor, voltage, target_torque)
    # Find slip required to maintain target torque at given voltage
    # Iteratively solve for slip that produces desired torque
    
    for s in slip_range
        p_m, p_e, q, i, torque, eff = calculate_motor_performance(motor, voltage, s)
        if abs(torque - target_torque) < tolerance
            return s
        end
    end
end

function analyze_constant_load_operation(motor, target_torque)
    voltage_range = range(0.6, 1.2, length=50)
    
    for voltage in voltage_range
        required_slip = calculate_motor_slip_for_load(motor, voltage, target_torque)
        r2_effective = motor.r2 / required_slip
        
        # Calculate performance at this operating point
        p_m, p_e, q, i, torque, efficiency = 
            calculate_motor_performance(motor, voltage, required_slip)
        
        # Observe: efficiency drops at low voltage due to higher slip
    end
end
```

## 8. Key Characteristics and Applications

### 8.1 U-P-Q Behavior

1. **U-shaped Q vs V curve**: Fundamental characteristic for power system studies
2. **Minimum reactive power**: Occurs at optimal voltage (≈0.8-1.0 p.u.)
3. **Power factor variation**: Changes significantly with voltage and load

### 8.2 Power System Impact

1. **Voltage regulation**: Motors affect system voltage through reactive power consumption
2