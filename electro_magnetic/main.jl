# main.jl - Main Electromagnetic Transient Simulation Controller
# AC Distribution System with Inverter-Based Resources (IBRs)

using DifferentialEquations, LinearAlgebra, Plots, Random, Statistics
using SparseArrays, NLsolve, FFTW, Polynomials

# Include all component files
include("parameters.jl")
include("components.jl")
include("transformations.jl")
include("models.jl")
include("control.jl")
include("system_dynamics.jl")
include("initialization.jl")
include("analysis.jl")
# include("visualization.jl")

# Set random seed for reproducibility
Random.seed!(42)

function main()
    println("=== Electromagnetic Transient Simulation ===")
    println("AC Distribution System with IBRs")
    
    # Initialize system components
    println("\nInitializing system components...")
    components = initialize_all_components()
    
    # Set initial conditions
    println("Setting initial conditions...")
    x0 = set_initial_conditions()
    
    # Setup simulation parameters
    tspan = (0.0, T_SIM)
    
    println("Starting electromagnetic transient simulation...")
    println("Simulation time: $(T_SIM) seconds")
    println("Time step: $(DT) seconds")
    
    # Create and solve ODE problem
    prob = ODEProblem(system_dynamics!, x0, tspan, components)
    
    # Use Rodas5 solver for stiff systems (electromagnetic transients)
    # Disable autodiff to avoid issues with non-compatible functions
    sol = solve(prob, Rodas5(autodiff=false), dt=DT, adaptive=true, reltol=1e-6, abstol=1e-8)
    
    println("Simulation completed successfully!")
    
    # Analyze results
    println("\nAnalyzing results...")
    results = analyze_simulation_results(sol)
    
    # Create visualizations
    println("Creating visualization plots...")
    create_all_plots(sol, results)
    
    # Generate report
    generate_simulation_report(results)
    
    return sol, results
end

# Run simulation if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    sol, results = main()
end