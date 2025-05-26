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
    
    # Thevenin equivalent circuit
    Zth = (Zm * Z1) / (Zm + Z1)
    Vth = V1 * Zm / (Zm + Z1)
    
    # Total impedance
    Ztotal = Zth + Z2
    
    # Rotor current (referred to stator)
    I2 = Vth / Ztotal
    
    # Stator current
    Im = V1 / Zm  # Magnetizing current
    I1 = I2 + Im
    
    # Power calculations (3-phase)
    P_airgap = 3 * abs(I2)^2 * motor.r2 / slip  # Air gap power
    P_rotor_loss = 3 * abs(I2)^2 * motor.r2     # Rotor copper loss
    P_mech = P_airgap - P_rotor_loss             # Mechanical power
    P_stator_loss = 3 * abs(I1)^2 * motor.r1    # Stator copper loss
    P_elec = P_airgap + P_stator_loss            # Electrical power input
    
    # Reactive power
    Q = 3 * imag(V1 * conj(I1))
    
    # Torque calculation
    synchronous_speed = 120 * motor.rated_frequency / motor.poles  # rpm
    ωs = 2π * synchronous_speed / 60  # rad/s
    torque = P_airgap / ωs
    
    # Efficiency
    efficiency = P_mech > 0 ? P_mech / P_elec * 100 : 0  # percentage
    
    return P_mech/1000, P_elec/1000, Q/1000, abs(I1), torque, efficiency
end

"""
    plot_voltage_characteristics(motor::InductionMotor; voltage_range=(0.5, 1.2))

Plot motor characteristics vs voltage at rated slip.
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
    
    # Create subplots
    p1 = plot(voltage_pu, P_mech, 
              label="Mechanical Power", 
              xlabel="Voltage (p.u.)", 
              ylabel="Power (kW)",
              title="Power vs Voltage",
              linewidth=2,
              color=:blue)
    plot!(voltage_pu, P_elec, 
          label="Electrical Power", 
          linewidth=2,
          color=:red,
          linestyle=:dash)
    
    p2 = plot(voltage_pu, Q_reactive,
              label="Reactive Power",
              xlabel="Voltage (p.u.)",
              ylabel="Reactive Power (kVAr)",
              title="Reactive Power vs Voltage",
              linewidth=2,
              color=:green)
    
    p3 = plot(voltage_pu, current,
              label="Stator Current",
              xlabel="Voltage (p.u.)",
              ylabel="Current (A)",
              title="Current vs Voltage",
              linewidth=2,
              color=:orange)
    
    p4 = plot(voltage_pu, efficiency,
              label="Efficiency",
              xlabel="Voltage (p.u.)",
              ylabel="Efficiency (%)",
              title="Efficiency vs Voltage",
              linewidth=2,
              color=:purple)
    
    # Combine plots
    combined_plot = plot(p1, p2, p3, p4, 
                        layout=(2,2), 
                        size=(800, 600),
                        plot_title="Induction Motor Voltage Characteristics")
    
    return combined_plot
end

"""
    plot_pq_characteristics(motor::InductionMotor; voltage_levels=[0.8, 0.9, 1.0, 1.1, 1.2])

Plot P-Q characteristics for different voltage levels.
"""
function plot_pq_characteristics(motor::InductionMotor; voltage_levels=[0.8, 0.9, 1.0, 1.1, 1.2])
    pq_plot = plot(xlabel="Active Power (kW)", 
                   ylabel="Reactive Power (kVAr)",
                   title="Induction Motor P-Q Characteristics",
                   legend=:topright,
                   size=(700, 500))
    
    colors = palette(:viridis, length(voltage_levels))
    
    for (i, v_pu) in enumerate(voltage_levels)
        voltage = v_pu * motor.rated_voltage
        slip_range = range(0.001, 0.3, length=100)
        
        P_vals = Float64[]
        Q_vals = Float64[]
        
        for s in slip_range
            p_m, p_e, q, i, t, eff = calculate_motor_performance(motor, voltage, s)
            push!(P_vals, p_m)
            push!(Q_vals, q)
        end
        
        plot!(P_vals, Q_vals, 
              label="V = $(v_pu) p.u.",
              linewidth=2,
              color=colors[i])
        
        # Mark rated operating point
        if v_pu == 1.0
            synchronous_speed = 120 * motor.rated_frequency / motor.poles
            rated_slip = (synchronous_speed - motor.rated_speed) / synchronous_speed
            p_rated, _, q_rated, _, _, _ = calculate_motor_performance(motor, voltage, rated_slip)
            scatter!([p_rated], [q_rated], 
                    markersize=8, 
                    color=:red, 
                    label="Rated Point",
                    markershape=:star)
        end
    end
    
    return pq_plot
end

"""
    plot_3d_upq_surface(motor::InductionMotor)

Create a 3D surface plot showing the relationship between U, P, and Q.
"""
function plot_3d_upq_surface(motor::InductionMotor)
    # Define ranges
    voltage_range = range(0.6, 1.3, length=30)
    slip_range = range(0.01, 0.3, length=30)
    
    # Pre-allocate arrays
    U_grid = Float64[]
    P_grid = Float64[]
    Q_grid = Float64[]
    
    # Calculate performance for each combination
    for v_pu in voltage_range
        for s in slip_range
            voltage = v_pu * motor.rated_voltage
            p_m, _, q, _, _, _ = calculate_motor_performance(motor, voltage, s)
            
            push!(U_grid, v_pu)
            push!(P_grid, p_m)
            push!(Q_grid, q)
        end
    end
    
    # Create 3D surface plot
    surface_plot = scatter(U_grid, P_grid, Q_grid,
                          xlabel="Voltage (p.u.)",
                          ylabel="Active Power (kW)",
                          zlabel="Reactive Power (kVAr)",
                          title="3D U-P-Q Relationship for Induction Motor",
                          markersize=2,
                          markercolor=:viridis,
                          camera=(45, 30),
                          size=(800, 600))
    
    return surface_plot
end

"""
    demonstrate_induction_motor_analysis()

Comprehensive demonstration of induction motor analysis capabilities.
"""
function demonstrate_induction_motor_analysis()
    println("="^60)
    println("INDUCTION MOTOR ANALYSIS DEMONSTRATION")
    println("="^60)
    
    # Create example motor
    motor = InductionMotor(
        rated_power=15.0,      # 15 kW
        rated_voltage=400.0,   # 400 V
        rated_frequency=50.0,  # 50 Hz
        rated_speed=1450.0,    # 1450 rpm
        poles=4,
        r1=0.4,               # Stator resistance
        x1=0.8,               # Stator reactance
        r2=0.25,              # Rotor resistance
        x2=0.8,               # Rotor reactance
        xm=25.0               # Magnetizing reactance
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
    
    # Voltage characteristics
    p1 = plot_voltage_characteristics(motor)
    savefig(p1, "motor_voltage_characteristics.png")
    println("   ✓ Voltage characteristics plot saved")
    
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
    
    return motor, p1, p2, p3
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
    motor, p1, p2, p3 = demonstrate_induction_motor_analysis()
    
    # Create comparison plot
    println("\n5. Creating Motor Comparison...")
    p4 = create_motor_comparison()
    savefig(p4, "motor_size_comparison.png")
    println("   ✓ Motor size comparison plot saved")
    
    # Display all plots
    display(p1)
    display(p2)
    display(p3)
    display(p4)
    
    println("\n" * "="^60)
    println("INDUCTION MOTOR ANALYSIS COMPLETED")
    println("All plots saved and displayed successfully!")
    println("="^60)
    
    return motor, [p1, p2, p3, p4]
end

# Export functions for use in other modules
export InductionMotor, calculate_motor_performance
export plot_voltage_characteristics, plot_pq_characteristics, plot_3d_upq_surface
export demonstrate_induction_motor_analysis, create_motor_comparison

# Run main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
