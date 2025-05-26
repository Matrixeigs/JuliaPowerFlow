# Static Generator Mathematical Models

This document provides a comprehensive mathematical formulation of static generator models used in power system analysis, particularly for cylindrical rotor synchronous machines.

## 1. Basic Power Equations

For any cylindrical rotor synchronous machine, the active and reactive power outputs are governed by:

$$
P = \frac{E_{q}U_{t}}{X_{d}}\sin \delta
$$

$$
Q = \frac{E_{q}U_{t}\cos\delta - U_{t}^{2}}{X_{d}}
$$

where:
- $P$ is the active power output (p.u.)
- $Q$ is the reactive power output (p.u.)
- $E_q$ is the internal voltage magnitude behind direct-axis reactance (p.u.)
- $U_t$ is the terminal voltage magnitude (p.u.)
- $X_d$ is the direct-axis synchronous reactance (p.u.)
- $\delta$ is the power angle (radians)

## 2. Generator Capability Constraints

### 2.1 Apparent Power Limit

The fundamental thermal constraint limits the apparent power:

$$
P^{2} + Q^{2} \leq S_{\max}^{2}
$$

where $S_{\max}$ is the maximum apparent power rating (p.u.).

### 2.2 Active Power Limits

Mechanical constraints impose active power limits:

$$
P_{\min} \leq P \leq P_{\max}
$$

Typically:
- $P_{\min} = 0$ for thermal units (or minimum stable generation)
- $P_{\max}$ is the rated active power output

### 2.3 Power Angle Limits

Stability considerations limit the power angle:

$$
\delta_{\min} \leq \delta \leq \delta_{\max}
$$

Common values:
- $\delta_{\max} = 60°$ for stable operation
- $\delta_{\min} = -60°$ for motor operation (if applicable)

### 2.4 Internal Voltage Limits

Field winding constraints limit the internal voltage:

$$
E_{q,\min} \leq E_q \leq E_{q,\max}
$$

where:
- $E_{q,\min}$ corresponds to minimum field excitation
- $E_{q,\max}$ corresponds to maximum field excitation

### 2.5 Reactive Power Limits

Direct reactive power constraints:

$$
Q_{\min} \leq Q \leq Q_{\max}
$$

These may be constant or functions of active power output.

## 3. Capability Region Models

### 3.1 Simple Box Model

The simplest representation uses constant limits:

```julia
function simple_box_constraints(P, Q, P_min, P_max, Q_min, Q_max)
    return (P_min ≤ P ≤ P_max) && (Q_min ≤ Q ≤ Q_max)
end
```

### 3.2 Cylindrical Rotor Model

Combines apparent power and active power constraints:

```julia
function cylindrical_rotor_constraints(P, Q, P_min, P_max, S_max)
    power_limits = (P_min ≤ P ≤ P_max)
    apparent_limit = (P² + Q² ≤ S_max²)
    return power_limits && apparent_limit
end
```

### 3.3 Full Synchronous Machine Model

Incorporates all physical constraints:

```julia
function full_generator_constraints(P, Q, U_t, params)
    # Extract parameters
    P_min, P_max = params.P_min, params.P_max
    Q_min, Q_max = params.Q_min, params.Q_max
    S_max = params.S_max
    δ_max = params.δ_max
    E_q_min, E_q_max = params.E_q_min, params.E_q_max
    X_d = params.X_d
    
    # Check basic limits
    if !(P_min ≤ P ≤ P_max) || !(Q_min ≤ Q ≤ Q_max)
        return false
    end
    
    # Check apparent power limit
    if P² + Q² > S_max²
        return false
    end
    
    # Check synchronous machine constraints
    if P > 0  # Only for generating mode
        # Calculate required E_q and δ
        δ_rad = δ_max * π / 180
        E_q_required = P * X_d / (U_t * sin(δ_rad))
        
        if E_q_required < E_q_min || E_q_required > E_q_max
            return false
        end
    end
    
    return true
end
```

## 4. Dynamic Q_min Calculation

For realistic generator modeling, $Q_{\min}$ is often a function of active power:

$$
Q_{\min}(P) = \max\left(Q_{\min,\text{constant}}, \frac{E_{q,\min}U_t\cos(\delta_{\max}) - U_t^2}{X_d}\right)
$$

This accounts for the underexcitation limit that varies with power output.

### 4.1 Implementation

```julia
function calculate_qmin_function(P, U_t, E_q_min, δ_max, X_d, Q_min_constant)
    if P ≤ 0
        return Q_min_constant
    end
    
    δ_rad = δ_max * π / 180
    Q_from_underexcitation = (E_q_min * U_t * cos(δ_rad) - U_t²) / X_d
    
    return max(Q_min_constant, Q_from_underexcitation)
end
```

## 5. Capability Curve Construction

### 5.1 P-Q Plane Representation

The generator capability region in the P-Q plane is bounded by:

1. **Circular arc**: $P^2 + Q^2 = S_{\max}^2$ (thermal limit)
2. **Horizontal lines**: $P = P_{\max}$, $P = P_{\min}$ (turbine limits)
3. **Curved boundaries**: From synchronous machine constraints
4. **Vertical/sloped lines**: $Q = Q_{\max}$, $Q = Q_{\min}(P)$

### 5.2 Parametric Representation

For plotting capability curves, use parametric equations:

```julia
function generate_capability_curve(params, U_t=1.0, n_points=100)
    # Thermal limit (circular)
    θ = range(0, 2π, length=n_points)
    P_thermal = params.S_max .* cos.(θ)
    Q_thermal = params.S_max .* sin.(θ)
    
    # Power angle limits
    δ_range = range(-params.δ_max, params.δ_max, length=n_points) .* π/180
    P_delta_max = (params.E_q_max * U_t / params.X_d) .* sin.(δ_range)
    Q_delta_max = (params.E_q_max * U_t * cos.(δ_range) .- U_t²) / params.X_d
    
    P_delta_min = (params.E_q_min * U_t / params.X_d) .* sin.(δ_range)
    Q_delta_min = (params.E_q_min * U_t * cos.(δ_range) .- U_t²) / params.X_d
    
    return (P_thermal, Q_thermal), (P_delta_max, Q_delta_max), (P_delta_min, Q_delta_min)
end
```

## 6. PV to PQ Bus Transition

During power flow iterations, PV buses may violate reactive power limits and transition to PQ buses:

### 6.1 Transition Criteria

A PV bus transitions to PQ when:
- $Q_{\text{calculated}} > Q_{\max}$ → Set $Q = Q_{\max}$
- $Q_{\text{calculated}} < Q_{\min}(P)$ → Set $Q = Q_{\min}(P)$

### 6.2 Implementation Strategy

```julia
function check_pv_bus_limits(bus, Q_calculated, P_scheduled)
    Q_max = bus.generator.Q_max
    Q_min = calculate_qmin_function(bus.generator, P_scheduled)
    
    if Q_calculated > Q_max
        return :transition_to_pq, Q_max
    elseif Q_calculated < Q_min
        return :transition_to_pq, Q_min
    else
        return :remain_pv, Q_calculated
    end
end
```

## 7. Modeling Considerations

### 7.1 Voltage Dependency

For detailed studies, consider voltage dependency:

$$
E_q = E_q^0 + K_e(V_{\text{ref}} - V_t)
$$

where $K_e$ is the excitation system gain.

### 7.2 Frequency Effects

For frequency studies, include governor response:

$$
P = P^0 + K_g(f_{\text{ref}} - f)
$$

where $K_g$ is the governor droop constant.

### 7.3 Temperature Effects

For thermal units, consider temperature derating:

$$
P_{\max}(T) = P_{\max,\text{rated}} \cdot \left(1 - k_T(T - T_{\text{rated}})\right)
$$

## 8. Advanced Features

### 8.1 Multi-Machine Modeling

For systems with multiple generators at one bus:

$$
P_{\text{bus}} = \sum_{i=1}^{n} P_i, \quad Q_{\text{bus}} = \sum_{i=1}^{n} Q_i
$$

Each generator maintains its individual constraints.

### 8.2 Renewable Integration

For variable renewable sources:

$$
P_{\text{available}}(t) = P_{\text{rated}} \cdot CF(t)
$$

where $CF(t)$ is the capacity factor.

### 8.3 Energy Storage

For storage systems, add energy constraints:

$$
\int_0^t P(\tau) d\tau \leq E_{\max}
$$

## 9. Implementation Notes

### 9.1 Numerical Considerations

- Use appropriate tolerances for constraint checking
- Handle boundary cases carefully
- Consider using smooth approximations for sharp constraints

### 9.2 Convergence Enhancement

- Implement adaptive Q limits during iterations
- Use continuation methods for difficult cases
- Consider warm-starting from previous solutions

### 9.3 Validation

- Compare results with manufacturer data
- Validate against field measurements
- Cross-check with commercial software

This comprehensive mathematical framework provides the foundation for accurate static generator modeling in power system analysis applications.