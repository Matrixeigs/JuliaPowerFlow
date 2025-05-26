# JuliaPowerFlow
julia based power flow implementation

# Julia Power Flow - Generator Capability Curves

This repository contains Julia code to visualize and analyze generator capability curves used in power flow analysis.

## Static Generation Feasibility Regions

The `static_generation.jl` file provides functions to plot and visualize different types of generator capability curves in the QP-plane (Q on x-axis, P on y-axis):

1. **Simple Box Model**: Rectangular PQ region with independent P and Q limits
2. **Cylindrical Rotor Model**: Based on apparent power constraint ($P^2 + Q^2 \leq S_{max}^2$) with P limits
3. **Model with Delta Limits**: Incorporates power angle constraints ($\delta_{min} \leq \delta \leq \delta_{max}$)
4. **Full Constraints Model**: Includes all constraints (P limits, Q limits, S limit, delta limits, and E_q limits)

## Usage Examples

```julia
# Import the module
include("static_generation.jl")

# Plot individual models
plot_simple_box_model(p_min=-1.0, p_max=1.0, q_min=-0.8, q_max=0.8)
plot_cylindrical_rotor_model(p_min=-1.0, p_max=1.0, s_max=1.2)
plot_with_delta_limits(p_min=-1.0, p_max=1.0, s_max=1.0, delta_max=60, eq=1.2, ut=1.0, xd=0.8)
plot_full_constraints(p_min=-1.0, p_max=1.0, q_min=-0.8, q_max=0.8, s_max=1.0, delta_max=60, 
                      eq_min=0.8, eq_max=1.2, ut=1.0, xd=0.8)

# Compare all models
compare_all_models()

# Run demonstration with examples
demonstrate_models()
```

## Parameter Descriptions

- `p_max`: Maximum active power output (pu)
- `p_min`: Minimum active power output (pu)
- `q_max`: Maximum reactive power output (pu)
- `q_min`: Minimum reactive power output (pu)
- `s_max`: Maximum apparent power (pu)
- `delta_max`: Maximum power angle (degrees, typically ≤ 60°)
- `eq_min`: Minimum internal voltage magnitude (pu)
- `eq_max`: Maximum internal voltage magnitude (pu)
- `ut`: Terminal voltage magnitude (pu)
- `xd`: Direct-axis synchronous reactance (pu)

## Mathematical Model

The capability curves implement the constraints for cylindrical rotor synchronous machines:

$$
P^{2} + Q^{2} \leq S_{\max}^{2} \\
P = \frac{E_{q}U_{t}}{X_{d}}\sin \delta \\
Q = \frac{E_{q}U_{t}\cos\delta - U_{t}^{2}}{X_{d}} \\
P_{\min} \leq P \leq P_{\max} \\
\delta_{\min} \leq \delta \leq \delta_{\max} \\
E_{q,\min} \leq E_q \leq E_{q,\max} \\
Q_{\min} \leq Q \leq Q_{\max}
$$

## Requirements

- Julia 1.0 or higher
- Plots.jl package
