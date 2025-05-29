include("static_generation.jl")
include("test_generator_D_shape.jl")  # Include the D-shape functions
include("reactive_power_examples.jl")  # Include reactive power examples

# Plot individual models
plot_simple_box_model(p_min=-1.0, p_max=0.9, q_min=-0.8, q_max=0.8)
plot_cylindrical_rotor_model(p_min=-1.0, p_max=0.9, s_max=1.2)
plot_with_delta_limits(p_min=-1.0, p_max=0.9, s_max=1.0, delta_max=60, eq=1.2, ut=1.0, xd=0.8)
plot_full_constraints(p_min=-1.0, p_max=0.9, q_min=-0.8, q_max=0.8, s_max=1.0, delta_max=60, 
                      eq_min=0.8, eq_max=1.2, ut=1.0, xd=0.8)

# Compare all models
compare_all_models()

# Run demonstration with examples
demonstrate_models()

# Test Q_min as function of P
println("Testing Q_min as function of P...")
qmin_func_plot = plot_qmin_function(
    p_min=0.0, p_max=0.9,
    eq_min=1.1, ut=1.0, xd=0.8, delta_max=45,
    title="Q_min as Function of P"
)
display(qmin_func_plot)

# Test PV→PQ transition scenarios
println("Testing PV→PQ transition scenarios...")
pv_pq_results = test_pv_pq_scenarios()
scenario1, scenario2, scenario3, combined = pv_pq_results

display(scenario1)
display(scenario2) 
display(scenario3)
display(combined)

# Test full constraints with Q_min function
println("Testing full constraints with Q_min function...")
full_with_qmin_func = plot_full_constraints(
    p_min=0.0, p_max=0.9, 
    q_min=-0.5, q_max=0.8, 
    delta_max=50,
    eq_min=1.1, eq_max=1.4, 
    ut=1.0, xd=0.8,
    use_qmin_function=true,
    title="Full Constraints with Q_min(P)"
)
display(full_with_qmin_func)

# Compare with constant Q_min
full_with_constant_qmin = plot_full_constraints(
    p_min=0.0, p_max=0.9, 
    q_min=-0.5, q_max=0.8, 
    s_max=1.0, 
    delta_max=50,
    eq_min=1.1, eq_max=1.4, 
    ut=1.0, xd=0.8,
    use_qmin_function=false,
    title="Full Constraints with Constant Q_min"
)
display(full_with_constant_qmin)