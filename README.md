# JuliaPowerFlow

A comprehensive power system analysis package written in Julia, featuring Newton-Raphson power flow analysis, generator capability modeling, and advanced visualization capabilities.

## 📁 Project Structure

```
JuliaPowerFlow/
├── src/                          # Source code
│   └── JuliaPowerFlow.jl        # Main module
├── data/                         # Test case data
│   └── case9.jl                 # IEEE 9-bus system
├── algorithms/                   # Core algorithms
│   └── power_flow.jl            # Newton-Raphson power flow
├── visualization/                # Visualization tools
│   └── power_flow_visualization.jl
├── examples/                     # Usage examples
│   ├── usage_example.jl         # Main usage demonstration
│   └── basic_newton_raphson.jl  # Newton-Raphson examples
├── tests/                        # Test suite
│   └── test_all_functions.jl    # Comprehensive tests
├── docs/                         # Documentation
│   └── static_generator.md      # Generator modeling docs
├── Project.toml                  # Package configuration
└── README.md                     # This file
```

## 🚀 Quick Start

1. **Load the package:**
```julia
push!(LOAD_PATH, "src")
using JuliaPowerFlow
```

2. **Run power flow analysis:**
```julia
case_data = case9()
V, S, Sf, St = run_power_flow(case_data)
```

3. **Create power system components:**
```julia
sys = PowerSystem(100.0)  # 100 MVA base
bus = Bus(1, bus_type=SLACK_BUS, voltage_magnitude=1.04)
gen = Generator(1, 1, p_max=250.0, s_max=300.0)
add_component!(sys, bus)
add_component!(sys, gen)
```

## 📊 Features

### Power Flow Analysis
- Newton-Raphson solver with automatic convergence monitoring
- Support for PV, PQ, and slack buses
- Line flow calculations and loss analysis
- Comprehensive result visualization

### Generator Modeling
- Cylindrical rotor synchronous machine constraints
- P-Q capability region analysis
- Dynamic Q_min calculation based on power angle limits
- Generator capability curve visualization

### Visualization
- Interactive power flow convergence plots
- Generator capability curve plotting
- System topology visualization
- Jacobian matrix heatmaps with animation

### Test Cases
- IEEE 9-bus standard test system
- Modular component-based system building
- Extensible framework for custom test cases

## 📖 Usage Examples

### Basic Power Flow
```julia
# Load IEEE 9-bus system
case_data = case9()

# Run power flow analysis
V, S, Sf, St = run_power_flow(case_data)

# Run with visualization
V, S, Sf, St = run_power_flow_with_visualization(case_data)
```

### Generator Capability Analysis
```julia
# Create generator
gen = Generator(1, 1, p_max=250.0, q_max=100.0, s_max=300.0)

# Test operating point
valid, msg = check_capability_constraints(gen, 200.0, 80.0)

# Calculate Q_min as function of P
q_min = calculate_qmin_function(gen, 150.0)
```

### Component-Based System Building
```julia
# Create system
sys = PowerSystem(100.0)

# Add components
add_component!(sys, Bus(1, bus_type=SLACK_BUS))
add_component!(sys, Generator(1, 1, p_max=100.0))
add_component!(sys, Branch(1, 1, 2, reactance=0.1))

# Get system information
println("Buses: $(get_bus_count(sys))")
```

## 🧪 Running Tests

```julia
# Run comprehensive test suite
include("tests/test_all_functions.jl")

# Run specific examples
include("examples/basic_newton_raphson.jl")
include("examples/usage_example.jl")
```

## 📚 Documentation

- **Generator Modeling**: See `docs/static_generator.md` for detailed mathematical formulation
- **API Reference**: All functions include comprehensive docstrings
- **Examples**: Check `examples/` directory for usage demonstrations

## 🛠️ Development

The package is organized into modular components:
- **Algorithms**: Core computational methods
- **Data**: Test case definitions and system data
- **Visualization**: Plotting and animation tools
- **Examples**: Demonstration scripts
- **Tests**: Validation and verification

## 🔧 Dependencies

- LinearAlgebra.jl
- SparseArrays.jl
- Plots.jl
- Printf.jl
- Statistics.jl
- ColorSchemes.jl

## 📄 License

This project is developed for educational and research purposes.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.
