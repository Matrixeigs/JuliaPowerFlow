# Synchronous Generator Seventh-Order Mathematical Model

## Overview

This document describes the mathematical model of a synchronous generator using the seventh-order representation. This model captures the electromagnetic transients in both d-axis and q-axis circuits, including transient and subtransient phenomena.

## State Variables

The model uses seven state variables:

| Variable | Symbol | Description | Unit |
|----------|--------|-------------|------|
| δ | δ | Rotor angle | rad |
| ω | ω | Rotor speed deviation from synchronous | pu |
| E'q | E'q | q-axis transient EMF | pu |
| E'd | E'd | d-axis transient EMF | pu |
| E''q | E''q | q-axis subtransient EMF | pu |
| E''d | E''d | d-axis subtransient EMF | pu |
| ψ2q | ψ2q | q-axis damper winding flux linkage | pu |

## Differential Equations

### 1. Rotor Dynamics (Swing Equations)

**Rotor angle equation:**
```
dδ/dt = ωs × ω
```

**Rotor speed equation:**
```
dω/dt = (Pm - Pe - D×ω) / (2H)
```

Where:
- ωs = synchronous speed (377 rad/s for 60 Hz)
- Pm = mechanical power input (pu)
- Pe = electrical power output (pu)
- D = damping coefficient (pu)
- H = inertia constant (s)

### 2. d-axis Transient EMF

```
dE'q/dt = (Efd - E'q + (Xd - X'd)×Id) / T'd0
```

Where:
- Efd = field voltage (pu)
- Xd = d-axis synchronous reactance (pu)
- X'd = d-axis transient reactance (pu)
- T'd0 = d-axis transient open-circuit time constant (s)
- Id = d-axis stator current (pu)

### 3. q-axis Transient EMF

```
dE'd/dt = (-E'd - (Xq - X'q)×Iq) / T'q0
```

Where:
- Xq = q-axis synchronous reactance (pu)
- X'q = q-axis transient reactance (pu)
- T'q0 = q-axis transient open-circuit time constant (s)
- Iq = q-axis stator current (pu)

### 4. d-axis Subtransient EMF

```
dE''q/dt = (E'q - E''q + (X'd - X''d)×Id) / T''d
```

Where:
- X''d = d-axis subtransient reactance (pu)
- T''d = d-axis subtransient short-circuit time constant (s)

### 5. q-axis Subtransient EMF

```
dE''d/dt = (E'd - E''d - (X'q - X''q)×Iq) / T''q
```

Where:
- X''q = q-axis subtransient reactance (pu)
- T''q = q-axis subtransient short-circuit time constant (s)

### 6. q-axis Damper Winding Flux

```
dψ2q/dt = (-ψ2q + E'd) / T''q0
```

Where:
- T''q0 = q-axis subtransient open-circuit time constant (s)

## Network Interface Equations

### Terminal Voltage Components

```
Vd = Vt × sin(δ)
Vq = Vt × cos(δ)
```

Where:
- Vt = terminal voltage magnitude (pu)

### Stator Current Calculation

The stator currents are calculated from the network equations:

```
Id = (E''q - Vq) / (X''d + Xe)
Iq = -(E''d - Vd) / (X''q + Xe)
```

Where:
- Xe = external reactance (transmission line reactance) (pu)

### Electrical Power Output

```
Pe = Vd × Id + Vq × Iq
```

### Electromagnetic Torque

```
Te = E''q × Iq + E''d × Id
```

## Parameter Definitions

### Machine Parameters

| Parameter | Symbol | Typical Value | Description |
|-----------|--------|---------------|-------------|
| H | H | 3-6 s | Inertia constant |
| D | D | 1-2 pu | Damping coefficient |
| Xd | Xd | 1.5-2.0 pu | d-axis synchronous reactance |
| Xq | Xq | 1.4-1.8 pu | q-axis synchronous reactance |
| X'd | X'd | 0.2-0.4 pu | d-axis transient reactance |
| X'q | X'q | 0.4-0.8 pu | q-axis transient reactance |
| X''d | X''d | 0.15-0.3 pu | d-axis subtransient reactance |
| X''q | X''q | 0.15-0.3 pu | q-axis subtransient reactance |

### Time Constants

| Parameter | Symbol | Typical Value | Description |
|-----------|--------|---------------|-------------|
| T'd0 | T'd0 | 5-10 s | d-axis transient open-circuit time constant |
| T'q0 | T'q0 | 0.5-2 s | q-axis transient open-circuit time constant |
| T''d0 | T''d0 | 0.02-0.05 s | d-axis subtransient open-circuit time constant |
| T''q0 | T''q0 | 0.02-0.1 s | q-axis subtransient open-circuit time constant |
| T''d | T''d | 0.02-0.05 s | d-axis subtransient short-circuit time constant |
| T''q | T''q | 0.02-0.1 s | q-axis subtransient short-circuit time constant |

## Initial Conditions Calculation

The initial conditions are calculated assuming steady-state operation:

### 1. Power Flow Solution

Given initial active power P0 and terminal voltage Vt:

```
δ0 = arcsin(P0 × X'd / (E'q0 × Vt))
```

### 2. Terminal Voltage Components

```
Vd0 = Vt × sin(δ0)
Vq0 = Vt × cos(δ0)
```

### 3. Initial Currents

From the power equation: S = P + jQ = Vt × I*

```
Id0 = Re(I × e^(j(π/2 - δ0)))
Iq0 = Im(I × e^(j(π/2 - δ0)))
```

### 4. Initial EMF Values

```
E'q0 = Vq0 + X'd × Id0
E'd0 = Vd0 - X'q × Iq0
E''q0 = Vq0 + X''d × Id0
E''d0 = Vd0 - X''q × Iq0
ψ2q0 = E'd0
```

## Model Assumptions

1. **Symmetrical three-phase system**: Only positive sequence components are considered
2. **Park's transformation**: dq0 reference frame aligned with the rotor
3. **Constant terminal voltage magnitude**: The infinite bus assumption
4. **Linear magnetic circuits**: Saturation effects are neglected
5. **Constant parameters**: Temperature and frequency dependencies are ignored
6. **No mechanical coupling**: Shaft dynamics are represented by simple swing equation

## Applications

This model is suitable for:

- Power system stability studies
- Transient analysis following faults or load changes
- Controller design (AVR, PSS)
- Small-signal stability analysis
- Electromagnetic transient studies (up to a few seconds)

## Limitations

- Not suitable for very fast electromagnetic transients (< 1 ms)
- Does not include magnetic saturation effects
- Assumes constant network impedance
- Does not model detailed transformer or transmission line dynamics
