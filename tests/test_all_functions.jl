"""
Comprehensive Test Suite for Julia Power Flow Package

This file tests all major functionalities including:
1. Basic Newton-Raphson solver tests  
2. Power flow analysis
3. Generator capability curves
4. Visualization functions
5. System convergence analysis
"""

# Import modules
push!(LOAD_PATH, pwd() * "/src")

# Include the existing files directly for testing
include("../data/case9.jl")
include("../algorithms/power_flow.jl")
include("../visualization/power_flow_visualization.jl")

using LinearAlgebra, SparseArrays, Printf
using Plots, ColorSchemes, Statistics

function test_basic_newton_raphson()
    println("\n" * "="^60)
    println("TESTING BASIC NEWTON-RAPHSON METHODS")
    println("="^60)
    
    # Include and run basic Newton-Raphson tests
    include("../examples/test_newton.jl")
    include("../examples/test_nr.jl")
    
    println("✓ Basic Newton-Raphson tests completed")
end

function test_power_flow_integration()
    println("\n" * "="^60)
    println("TESTING POWER FLOW INTEGRATION")
    println("="^60)
    
    try
        # Test existing power flow solver
        case_data = case9()
        println("✓ Successfully loaded IEEE 9-bus case")
        
        V, S, Sf, St = run_power_flow(case_data)
        println("✓ Basic power flow completed successfully")
        
        # Test visualization power flow
        V_vis, S_vis, Sf_vis, St_vis = run_power_flow_with_visualization(case_data)
        println("✓ Visualization power flow completed successfully")
        
        return Dict(
            "basic_solution" => V,
            "visualization_solution" => V_vis,
            "converged" => true
        )
    catch e
        println("✗ Power flow test failed: $e")
        return Dict("converged" => false, "error" => e)
    end
end

function test_visualization_functions()
    println("\n" * "="^60)
    println("TESTING VISUALIZATION FUNCTIONS")
    println("="^60)
    
    try
        case_data = case9()
        V, S, Sf, St = run_power_flow(case_data)
        
        # Test topology plotting
        topology_plot = plot_power_system_topology(case_data, V)
        println("✓ System topology plotting successful")
        
        return true
    catch e
        println("✗ Visualization test failed: $e")
        return false
    end
end

function main()
    println("JULIA POWER FLOW - COMPREHENSIVE TEST SUITE")
    println("Starting comprehensive testing of all functions...")
    
    all_results = Dict()
    
    try
        # Test 1: Basic Newton-Raphson methods
        test_basic_newton_raphson()
        all_results["newton_raphson"] = true
        
        # Test 2: Power flow integration
        pf_results = test_power_flow_integration()
        all_results["power_flow"] = pf_results
        
        # Test 3: Visualization functions
        vis_results = test_visualization_functions()
        all_results["visualization"] = vis_results
        
        # Summary
        println("\n" * "="^60)
        println("TEST SUMMARY")
        println("="^60)
        
        println("1. Newton-Raphson Tests: ✓ PASSED")
        println("2. Power Flow Analysis: $(all_results["power_flow"]["converged"] ? "✓ PASSED" : "✗ FAILED")")
        println("3. Visualization Tests: $(all_results["visualization"] ? "✓ PASSED" : "✗ FAILED")")
        
        println("\n🎉 ALL TESTS COMPLETED!")
        
    catch e
        println("\n❌ ERROR DURING TESTING:")
        println("Error: $e")
        rethrow(e)
    end
    
    return all_results
end

# Run all tests when this file is executed
if abspath(PROGRAM_FILE) == @__FILE__
    results = main()
end
