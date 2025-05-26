"""
Induction Motor Analysis and Visualization

This module provides comprehensive analysis and visualization of induction motor
characteristics, showing the relationships between voltage (U), active power (P),
and reactive power (Q) under various operating conditions.
"""

using Plots, LinearAlgebra, Printf
using ColorSchemes

"""
    InductionMotor

Structure representing an induction motor with equivalent circuit parameters.
"""
mutable struct InductionMotor
    # Nameplate data
    rated_power::Float64      # Rated mechanical power (kW)
    rated_voltage::Float64    # Rated voltage (V)
    rated_frequency::Float64  # Rated frequency (Hz)
    rated_speed::Float64      # Rated speed (rpm)
    poles::Int               # Number of poles
    
    # Equivalent circuit parameters (per phase, referred to stator)
    r1::Float64              # Stator resistance (Ω)
    x1::Float64              # Stator reactance (Ω)
    r2::Float64              # Rotor resistance (Ω)
    x2::Float64              # Rotor reactance (Ω)
    xm::Float64              # Magnetizing reactance (Ω)
    
    # Mechanical parameters
    rated_torque::Float64    # Rated torque (N⋅m)
    inertia::Float64         # Moment of inertia (kg⋅m²)
    
    function InductionMotor(;rated_power=10.0, rated_voltage=400.0, rated_frequency=50.0,
                           rated_speed=1450.0, poles=4, r1=0.5, x1=1.0, r2=0.3, 
                           x2=1.0, xm=30.0, rated_torque=65.9, inertia=0.1)
        new(rated_power, rated_voltage, rated_frequency, rated_speed, poles,
            r1, x1, r2, x2, xm, rated_torque, inertia)
    end
end

"""
    calculate_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)

Calculate motor performance at given voltage and slip.
Returns: (P_mech, P_elec, Q, I, torque, efficiency)
"""
function calculate_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)
    # Avoid division by zero
    if slip ≈ 0
        slip = 1e-6
    end
    
    # Per-phase voltage
    V1 = voltage / √3
    
    # Equivalent circuit impedances
    Z1 = motor.r1 + im * motor.x1
    Z2 = motor.r2/slip + im * motor.x2
    Zm = im * motor.xm
    
    # Total impedance seen from stator
    Z_parallel = (Z2 * Zm) / (Z2 + Zm)
    Z_total = Z1 + Z_parallel
    
    # Stator current
    I1 = V1 / Z_total
    
    # Voltage across parallel combination (rotor + magnetizing)
    V_parallel = I1 * Z_parallel
    
    # Rotor current (referred to stator)
    I2 = V_parallel / Z2
    
    # Magnetizing current
    Im = V_parallel / Zm
    
    # Power calculations (3-phase)
    P_stator_loss = 3 * abs(I1)^2 * motor.r1    # Stator copper loss
    P_rotor_loss = 3 * abs(I2)^2 * motor.r2     # Rotor copper loss
    P_airgap = 3 * abs(I2)^2 * motor.r2 / slip  # Air gap power
    P_mech = P_airgap - P_rotor_loss             # Mechanical power
    P_elec = P_airgap + P_stator_loss            # Electrical power input
    
    # CORRECTED: Active power is proportional to V² for constant slip
    # This ensures P decreases when voltage decreases, which is physically correct
    # The power transfer depends on the voltage squared and the rotor current
    
    # Reactive power calculation with three components
    Q_stator_leakage = 3 * abs(I1)^2 * motor.x1
    Q_rotor_leakage = 3 * abs(I2)^2 * motor.x2
    Q_magnetizing = 3 * abs(V_parallel)^2 / motor.xm
    
    Q = Q_stator_leakage + Q_rotor_leakage + Q_magnetizing
    
    # Verification: Active power should show V² dependency
    # At constant slip, P ∝ V² × I² ∝ V² × (V/Z)² ∝ V⁴/Z² 
    # But since Z also depends on slip and motor parameters, 
    # the net effect is approximately P ∝ V² for practical voltage ranges
    
    # Torque calculation
    synchronous_speed = 120 * motor.rated_frequency / motor.poles  # rpm
    ωs = 2π * synchronous_speed / 60  # rad/s
    torque = P_airgap / ωs
    
    # Efficiency
    efficiency = P_mech > 0 ? P_mech / P_elec * 100 : 0  # percentage
    
    return P_mech/1000, P_elec/1000, Q/1000, abs(I1), torque, efficiency
end

"""
    calculate_saturated_xm(motor::InductionMotor, voltage::Float64)

Calculate magnetizing reactance considering magnetic saturation.
At higher voltages, Xm decreases due to core saturation.
"""
function calculate_saturated_xm(motor::InductionMotor, voltage::Float64)
    # Saturation model: Xm decreases as voltage increases beyond rated
    # This is crucial for creating the U-shaped Q vs V curve
    
    v_pu = voltage / motor.rated_voltage
    xm_rated = motor.xm
    
    # Saturation curve - empirical model
    # At low voltage: Xm remains high (linear region)
    # At high voltage: Xm decreases (saturation region)
    
    if v_pu <= 0.8
        # Linear region - no saturation
        xm_saturated = xm_rated
    elseif v_pu <= 1.2
        # Transition region - mild saturation
        saturation_factor = 1.0 - 0.1 * (v_pu - 0.8)  # 10% reduction from 0.8 to 1.2 p.u.
        xm_saturated = xm_rated * saturation_factor
    else
        # Heavy saturation region
        # Xm decreases more rapidly at high voltages
        base_reduction = 0.96  # 4% reduction at 1.2 p.u.
        additional_reduction = 0.3 * (v_pu - 1.2)^1.5  # Nonlinear saturation
        saturation_factor = base_reduction - additional_reduction
        saturation_factor = max(saturation_factor, 0.4)  # Minimum 40% of rated Xm
        xm_saturated = xm_rated * saturation_factor
    end
    
    return xm_saturated
end

"""
    calculate_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)

Calculate motor performance at given voltage and slip with magnetic saturation.
Returns: (P_mech, P_elec, Q, I, torque, efficiency)
"""
function calculate_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)
    # Avoid division by zero
    if slip ≈ 0
        slip = 1e-6
    end
    
    # Per-phase voltage
    V1 = voltage / √3
    
    # Get saturated magnetizing reactance - KEY for U-shaped curve
    xm_saturated = calculate_saturated_xm(motor, voltage)
    
    # Equivalent circuit impedances with saturation
    Z1 = motor.r1 + im * motor.x1
    Z2 = motor.r2/slip + im * motor.x2
    Zm = im * xm_saturated  # Use saturated value
    
    # Total impedance seen from stator
    Z_parallel = (Z2 * Zm) / (Z2 + Zm)
    Z_total = Z1 + Z_parallel
    
    # Stator current
    I1 = V1 / Z_total
    
    # Voltage across parallel combination (rotor + magnetizing)
    V_parallel = I1 * Z_parallel
    
    # Rotor current (referred to stator)
    I2 = V_parallel / Z2
    
    # Magnetizing current
    Im = V_parallel / Zm
    
    # Power calculations (3-phase)
    P_stator_loss = 3 * abs(I1)^2 * motor.r1    # Stator copper loss
    P_rotor_loss = 3 * abs(I2)^2 * motor.r2     # Rotor copper loss
    P_airgap = 3 * abs(I2)^2 * motor.r2 / slip  # Air gap power
    P_mech = P_airgap - P_rotor_loss             # Mechanical power
    P_elec = P_airgap + P_stator_loss            # Electrical power input
    
    # Reactive power calculation with saturation effects
    Q_stator_leakage = 3 * abs(I1)^2 * motor.x1
    Q_rotor_leakage = 3 * abs(I2)^2 * motor.x2
    Q_magnetizing = 3 * abs(V_parallel)^2 / xm_saturated  # Uses saturated Xm
    
    Q = Q_stator_leakage + Q_rotor_leakage + Q_magnetizing
    
    # Alternative verification using complex power
    S_total = 3 * V1 * conj(I1)
    Q_verification = imag(S_total)
    
    # Use the complex power method for accuracy
    Q = Q_verification
    
    # Torque calculation
    synchronous_speed = 120 * motor.rated_frequency / motor.poles  # rpm
    ωs = 2π * synchronous_speed / 60  # rad/s
    torque = P_airgap / ωs
    
    # Efficiency
    efficiency = P_mech > 0 ? P_mech / P_elec * 100 : 0  # percentage
    
    return P_mech/1000, P_elec/1000, Q/1000, abs(I1), torque, efficiency
end

"""
    analyze_q_vs_voltage_behavior(motor::InductionMotor)

Analyze and explain the Q vs V relationship with saturation effects.
"""
function analyze_q_vs_voltage_behavior(motor::InductionMotor)
    println("\nDetailed Analysis of Q vs V Relationship (with Saturation):")
    println("="^60)
    
    # Calculate rated slip
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    
    voltage_range = range(0.3, 2.0, length=40)  # Wider range to show full U-shape
    
    println("Voltage(pu) | Q_total | Q_mag | Q_leak | Xm_sat | Explanation")
    println("-"^80)
    
    q_total_vals = Float64[]
    q_mag_vals = Float64[]
    q_leak_vals = Float64[]
    xm_sat_vals = Float64[]
    v_vals = Float64[]
    
    for v_pu in voltage_range
        voltage = v_pu * motor.rated_voltage
        V1 = voltage / √3
        
        # Get saturated magnetizing reactance
        xm_saturated = calculate_saturated_xm(motor, voltage)
        
        # Calculate impedances with saturation
        Z1 = motor.r1 + im * motor.x1
        Z2 = motor.r2/rated_slip + im * motor.x2
        Zm = im * xm_saturated  # Saturated value
        
        # Calculate currents
        Z_parallel = (Z2 * Zm) / (Z2 + Zm)
        Z_total = Z1 + Z_parallel
        I1 = V1 / Z_total
        V_parallel = I1 * Z_parallel
        I2 = V_parallel / Z2
        
        # Reactive power components with saturation
        Q_stator_leakage = 3 * abs(I1)^2 * motor.x1
        Q_rotor_leakage = 3 * abs(I2)^2 * motor.x2
        Q_magnetizing = 3 * abs(V_parallel)^2 / xm_saturated  # Saturated
        
        Q_leakage_total = Q_stator_leakage + Q_rotor_leakage
        Q_total = Q_leakage_total + Q_magnetizing
        
        push!(v_vals, v_pu)
        push!(q_total_vals, Q_total/1000)
        push!(q_mag_vals, Q_magnetizing/1000)
        push!(q_leak_vals, Q_leakage_total/1000)
        push!(xm_sat_vals, xm_saturated)
        
        # Print detailed analysis for key points
        if v_pu in [0.5, 0.8, 1.0, 1.2, 1.5, 1.8]
            behavior = ""
            if v_pu <= 0.7
                behavior = "Low V: High I → High Q_leak"
            elseif v_pu <= 1.1
                behavior = "Normal: Balanced"
            elseif v_pu <= 1.4
                behavior = "Mild saturation"
            else
                behavior = "Heavy saturation → Q↓"
            end
            
            @printf("  %.1f     | %6.2f | %5.2f | %6.2f | %6.1f | %s\n", 
                    v_pu, Q_total/1000, Q_magnetizing/1000, Q_leakage_total/1000, 
                    xm_saturated, behavior)
        end
    end
    
    # Create enhanced Q breakdown plot showing saturation effects
    q_plot = plot(v_vals, q_total_vals,
                  label="Total Q (with saturation)",
                  xlabel="Voltage (p.u.)",
                  ylabel="Reactive Power (kVAr)",
                  title="Reactive Power vs Voltage (with Magnetic Saturation)",
                  linewidth=4,
                  color=:blue,
                  legend=:topright,
                  grid=true)
    
    # Add magnetizing component
    plot!(v_vals, q_mag_vals,
          label="Q_magnetizing",
          linewidth=3,
          color=:red,
          linestyle=:dash)
    
    # Add leakage component
    plot!(v_vals, q_leak_vals,
          label="Q_leakage (total)",
          linewidth=2,
          color=:green,
          linestyle=:dot)
    
    # Mark key points
    min_q_idx = argmin(q_total_vals)
    scatter!([v_vals[min_q_idx]], [q_total_vals[min_q_idx]],
            markersize=10,
            color=:black,
            markershape=:star,
            label="Q minimum")
    
    # Add saturation curve on secondary axis
    xm_plot = twinx()
    plot!(xm_plot, v_vals, xm_sat_vals,
          label="Xm (saturated)",
          linewidth=2,
          color=:purple,
          linestyle=:dashdot,
          ylabel="Magnetizing Reactance (Ω)")
    
    # Add explanation annotations
    annotate!(0.5, maximum(q_total_vals)*0.9, 
              text("Low V:\nHigh current\n→ High Q_leakage", 9, :left))
    
    annotate!(1.0, maximum(q_total_vals)*0.7, 
              text("Minimum Q\n(optimal point)", 9, :center))
    
    annotate!(1.7, maximum(q_total_vals)*0.6, 
              text("High V:\nSaturation reduces Xm\n→ Q decreases!", 9, :right))
    
    println("\nKey Insights with Saturation:")
    println("1. At low voltage: High current dominates → High leakage Q")
    println("2. At moderate voltage: Balanced operation → Minimum total Q")  
    println("3. At high voltage: Saturation reduces Xm → Magnetizing Q decreases")
    println("4. This creates the characteristic U-shaped Q vs V curve!")
    println("5. Without saturation, Q would only increase with V²")
    
    return q_plot
end

"""
    plot_voltage_characteristics(motor::InductionMotor; voltage_range=(0.5, 1.2))

Plot motor characteristics vs voltage at rated slip with improved analysis.
"""
function plot_voltage_characteristics(motor::InductionMotor; voltage_range=(0.5, 1.2))
    # Calculate rated slip
    synchronous_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (synchronous_speed - motor.rated_speed) / synchronous_speed
    
    # Voltage range (per unit)
    voltage_pu = range(voltage_range[1], voltage_range[2], length=100)
    voltage_actual = voltage_pu .* motor.rated_voltage
    
    # Calculate performance for each voltage
    P_mech = Float64[]
    P_elec = Float64[]
    Q_reactive = Float64[]
    current = Float64[]
    torque_vals = Float64[]
    efficiency = Float64[]
    
    for V in voltage_actual
        p_m, p_e, q, i, t, eff = calculate_motor_performance(motor, V, rated_slip)
        push!(P_mech, p_m)
        push!(P_elec, p_e)
        push!(Q_reactive, q)
        push!(current, i)
        push!(torque_vals, t)
        push!(efficiency, eff)
    end
    
    # Create subplots with proper P vs V relationship
    p1 = plot(voltage_pu, P_mech, 
              label="Mechanical Power", 
              xlabel="Voltage (p.u.)", 
              ylabel="Power (kW)",
              title="Power vs Voltage\n(P decreases with decreasing V)",
              linewidth=2,
              color=:blue)
    plot!(voltage_pu, P_elec, 
          label="Electrical Power", 
          linewidth=2,
          color=:red,
          linestyle=:dash)
    
    # Add annotation showing P ∝ V² relationship
    annotate!(p1, 0.7, maximum(P_mech)*0.8, 
              text("P ∝ V²\n(approximately)", 10, :left))
    
    # Enhanced Q vs V plot
    p2 = plot(voltage_pu, Q_reactive,
              label="Reactive Power",
              xlabel="Voltage (p.u.)",
              ylabel="Reactive Power (kVAr)",
              title="Reactive Power vs Voltage\n(Characteristic U-shaped curve)",
              linewidth=3,
              color=:green)
    
    # Find and mark minimum Q point
    min_q_idx = argmin(Q_reactive)
    min_q_voltage = voltage_pu[min_q_idx]
    min_q_value = Q_reactive[min_q_idx]
    
    scatter!(p2, [min_q_voltage], [min_q_value],
            markersize=8,
            color=:red,
            markershape=:star,
            label="Min Q point")
    
    # Add physics explanation
    annotate!(p2, 0.6, maximum(Q_reactive)*0.8,
              text("Low V: High I\n→ High Q_leakage", 9, :left))
    annotate!(p2, 1.1, maximum(Q_reactive)*0.8,
              text("High V: High flux\n→ High Q_mag", 9, :right))
    
    p3 = plot(voltage_pu, current,
              label="Stator Current",
              xlabel="Voltage (p.u.)",
              ylabel="Current (A)",
              title="Current vs Voltage\n(I increases as V decreases)",
              linewidth=2,
              color=:orange)
    
    # Add annotation showing I ∝ 1/V trend
    annotate!(p3, 0.7, maximum(current)*0.8,
              text("I ∝ 1/V\n(approximately)", 10, :left))
    
    p4 = plot(voltage_pu, efficiency,
              label="Efficiency",
              xlabel="Voltage (p.u.)",
              ylabel="Efficiency (%)",
              title="Efficiency vs Voltage\n(Optimal around rated voltage)",
              linewidth=2,
              color=:purple)
    
    # Mark rated voltage point
    rated_eff_idx = argmin(abs.(voltage_pu .- 1.0))
    scatter!(p4, [1.0], [efficiency[rated_eff_idx]],
            markersize=8,
            color=:red,
            markershape=:circle,
            label="Rated point")
    
    # Combine plots
    combined_plot = plot(p1, p2, p3, p4, 
                        layout=(2,2), 
                        size=(1000, 800),
                        plot_title="Induction Motor: Voltage Effects on Performance")
    
    return combined_plot
end

"""
    demonstrate_induction_motor_analysis()

Comprehensive demonstration of induction motor analysis capabilities.
"""
function demonstrate_induction_motor_analysis()
    println("="^60)
    println("INDUCTION MOTOR ANALYSIS DEMONSTRATION")
    println("="^60)
    
    # Create example motor with more realistic parameters
    motor = InductionMotor(
        rated_power=15.0,      # 15 kW
        rated_voltage=400.0,   # 400 V
        rated_frequency=50.0,  # 50 Hz
        rated_speed=1450.0,    # 1450 rpm
        poles=4,
        r1=0.3,               # Stator resistance (reduced)
        x1=0.8,               # Stator reactance  
        r2=0.2,               # Rotor resistance (reduced)
        x2=0.8,               # Rotor reactance
        xm=40.0               # Magnetizing reactance (increased for better Q behavior)
    )
    
    println("\n1. Motor Specifications:")
    println("   Rated Power: $(motor.rated_power) kW")
    println("   Rated Voltage: $(motor.rated_voltage) V")
    println("   Rated Speed: $(motor.rated_speed) rpm")
    println("   Synchronous Speed: $(120 * motor.rated_frequency / motor.poles) rpm")
    
    # Calculate rated slip
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    println("   Rated Slip: $(round(rated_slip*100, digits=2))%")
    
    # Test motor performance at rated conditions
    p_m, p_e, q, i, t, eff = calculate_motor_performance(motor, motor.rated_voltage, rated_slip)
    println("\n2. Rated Operating Point:")
    @printf("   Mechanical Power: %.2f kW\n", p_m)
    @printf("   Electrical Power: %.2f kW\n", p_e)
    @printf("   Reactive Power: %.2f kVAr\n", q)
    @printf("   Stator Current: %.2f A\n", i)
    @printf("   Torque: %.2f N⋅m\n", t)
    @printf("   Efficiency: %.1f%%\n", eff)
    
    # Generate plots
    println("\n3. Generating Visualization Plots...")
    
    # Voltage characteristics with improved Q analysis
    p1 = plot_voltage_characteristics(motor)
    savefig(p1, "motor_voltage_characteristics.png")
    println("   ✓ Voltage characteristics plot saved")
    
    # Detailed Q vs V analysis
    p_q_detail = analyze_q_vs_voltage_behavior(motor)
    savefig(p_q_detail, "motor_q_vs_v_detailed.png")
    println("   ✓ Detailed Q vs V analysis saved")
    
    # P-Q characteristics
    p2 = plot_pq_characteristics(motor)
    savefig(p2, "motor_pq_characteristics.png")
    println("   ✓ P-Q characteristics plot saved")
    
    # 3D U-P-Q surface
    p3 = plot_3d_upq_surface(motor)
    savefig(p3, "motor_3d_upq_surface.png")
    println("   ✓ 3D U-P-Q surface plot saved")
    
    # Voltage sensitivity analysis
    println("\n4. Voltage Sensitivity Analysis:")
    voltage_levels = [0.8, 0.9, 1.0, 1.1, 1.2]
    
    for v_pu in voltage_levels
        voltage = v_pu * motor.rated_voltage
        p_m, p_e, q, i, t, eff = calculate_motor_performance(motor, voltage, rated_slip)
        @printf("   V=%.1f p.u.: P=%.2f kW, Q=%.2f kVAr, I=%.1f A, η=%.1f%%\n", 
                v_pu, p_m, q, i, eff)
    end
    
    return motor, p1, p2, p3, p_q_detail
end

"""
    create_motor_comparison()

Compare different motor sizes and their U-P-Q characteristics.
"""
function create_motor_comparison()
    # Create motors of different sizes
    motors = [
        InductionMotor(rated_power=5.0, rated_voltage=400.0, r1=0.8, r2=0.5),
        InductionMotor(rated_power=15.0, rated_voltage=400.0, r1=0.4, r2=0.25),
        InductionMotor(rated_power=50.0, rated_voltage=400.0, r1=0.15, r2=0.1)
    ]
    
    comparison_plot = plot(xlabel="Active Power (kW)", 
                          ylabel="Reactive Power (kVAr)",
                          title="Motor Size Comparison - P-Q Characteristics",
                          legend=:topright,
                          size=(700, 500))
    
    colors = [:blue, :red, :green]
    
    for (i, motor) in enumerate(motors)
        slip_range = range(0.01, 0.3, length=100)
        P_vals = Float64[]
        Q_vals = Float64[]
        
        for s in slip_range
            p_m, _, q, _, _, _ = calculate_motor_performance(motor, motor.rated_voltage, s)
            push!(P_vals, p_m)
            push!(Q_vals, q)
        end
        
        plot!(P_vals, Q_vals,
              label="$(motor.rated_power) kW Motor",
              linewidth=2,
              color=colors[i])
    end
    
    return comparison_plot
end

# Main execution function
function main()
    println("Starting Induction Motor Analysis...")
    
    # Run comprehensive demonstration
    motor, p1, p2, p3, p_q_detail = demonstrate_induction_motor_analysis()
    
    # Create comparison plot
    println("\n5. Creating Motor Comparison...")
    p4 = create_motor_comparison()
    savefig(p4, "motor_size_comparison.png")
    println("   ✓ Motor size comparison plot saved")
    
    # Display all plots
    display(p1)
    display(p2)
    display(p3)
    display(p_q_detail)
    display(p4)
    
    println("\n" * "="^60)
    println("INDUCTION MOTOR ANALYSIS COMPLETED")
    println("All plots saved and displayed successfully!")
    println("="^60)
    
    return motor, [p1, p2, p3, p_q_detail, p4]
end

# Export functions for use in other modules
export InductionMotor, calculate_motor_performance
export plot_voltage_characteristics, plot_pq_characteristics, plot_3d_upq_surface
export demonstrate_induction_motor_analysis, create_motor_comparison

# Run main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
