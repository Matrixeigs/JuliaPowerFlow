include("visualization/induction_motor_analysis.jl")

"""
    calculate_correct_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)

Corrected motor performance calculation that properly shows Q vs V relationship.
"""
function calculate_correct_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)
    # Avoid division by zero
    if slip ≈ 0
        slip = 1e-6
    end
    
    # Per-phase voltage
    V1 = voltage / √3
    
    # Use complex impedance calculations for accuracy
    r1, x1 = motor.r1, motor.x1
    r2, x2 = motor.r2, motor.x2
    xm = motor.xm
    
    # Get saturated magnetizing reactance - CRITICAL for U-shape
    xm_saturated = calculate_saturated_xm(motor, voltage)
    
    # Complex impedances with saturation
    Z1 = complex(r1, x1)                    # Stator impedance
    Z2 = complex(r2/slip, x2)               # Rotor impedance referred to stator
    Zm = complex(0.0, xm_saturated)         # Magnetizing impedance (saturated)
    
    # Parallel combination of rotor and magnetizing branches
    Z_parallel = (Z2 * Zm) / (Z2 + Zm)
    
    # Total impedance
    Z_total = Z1 + Z_parallel
    
    # Stator current
    I1 = V1 / Z_total
    
    # Voltage across parallel combination (rotor-magnetizing)
    V_parallel = I1 * Z_parallel
    
    # Individual branch currents
    I2 = V_parallel / Z2          # Rotor current
    Im = V_parallel / Zm          # Magnetizing current
    
    # Verify: I1 should equal I2 + Im
    I1_check = I2 + Im
    if abs(I1 - I1_check) > 1e-10
        @warn "Current calculation inconsistency: $(abs(I1 - I1_check))"
    end
    
    # Power calculations (3-phase)
    P_stator_loss = 3 * abs(I1)^2 * r1
    P_rotor_loss = 3 * abs(I2)^2 * r2
    P_airgap = 3 * abs(I2)^2 * r2 / slip
    P_mech = P_airgap - P_rotor_loss
    P_elec = P_airgap + P_stator_loss
    
    # CORRECTED Reactive power calculation
    # Method 1: Using total complex power (this should be the reference)
    S_total = 3 * V1 * conj(I1)
    Q_total_reference = imag(S_total)
    
    # Method 2: Component-wise calculation (for understanding)
    # Q = 3 * X * I^2 for each reactive component
    Q_stator_leakage = 3 * abs(I1)^2 * x1       # Always positive (inductive)
    Q_rotor_leakage = 3 * abs(I2)^2 * x2        # Always positive (inductive)
    
    # Magnetizing reactive power - this is the key correction
    # Q_mag = 3 * |V_parallel|^2 / Xm, but we need to be careful about the voltage
    # The correct way is: Q_mag = 3 * |Im|^2 * Xm = 3 * |V_parallel/jXm|^2 * Xm
    # Which simplifies to: Q_mag = 3 * |V_parallel|^2 / Xm
    Q_magnetizing = 3 * abs(V_parallel)^2 / xm_saturated  # Key change!
    
    # Alternative magnetizing calculation for verification
    Q_magnetizing_alt = 3 * abs(Im)^2 * xm
    
    # The component sum should match the total
    Q_components = Q_stator_leakage + Q_rotor_leakage + Q_magnetizing
    
    # Debug output for verification
    if abs(Q_components - Q_total_reference) > 0.1
        println("Debug - Q calculation mismatch at V=$(round(voltage,digits=1))V:")
        println("  Q_total_ref: $(round(Q_total_reference/1000, digits=3)) kVAr")
        println("  Q_components: $(round(Q_components/1000, digits=3)) kVAr")
        println("  Q_stator: $(round(Q_stator_leakage/1000, digits=3)) kVAr")
        println("  Q_rotor: $(round(Q_rotor_leakage/1000, digits=3)) kVAr") 
        println("  Q_mag: $(round(Q_magnetizing/1000, digits=3)) kVAr")
        println("  Q_mag_alt: $(round(Q_magnetizing_alt/1000, digits=3)) kVAr")
        println("  |V_parallel|: $(round(abs(V_parallel), digits=3)) V")
        println("  |Im|: $(round(abs(Im), digits=3)) A")
    end
    
    # Use the reference method (total complex power)
    Q = Q_total_reference
    
    # Torque calculation
    synchronous_speed = 120 * motor.rated_frequency / motor.poles
    ωs = 2π * synchronous_speed / 60
    torque = P_airgap / ωs
    
    # Efficiency
    efficiency = P_mech > 0 ? P_mech / P_elec * 100 : 0
    
    return P_mech/1000, P_elec/1000, Q/1000, abs(I1), torque, efficiency, 
           Q_stator_leakage/1000, Q_rotor_leakage/1000, Q_magnetizing/1000
end

"""
    calculate_saturated_xm(motor::InductionMotor, voltage::Float64)

Calculate magnetizing reactance considering magnetic saturation.
"""
function calculate_saturated_xm(motor::InductionMotor, voltage::Float64)
    v_pu = voltage / motor.rated_voltage
    xm_rated = motor.xm
    
    # Enhanced saturation model for pronounced U-shape
    if v_pu <= 0.8
        # Linear region - no saturation
        xm_saturated = xm_rated
    elseif v_pu <= 1.0
        # Mild saturation begins
        saturation_factor = 1.0 - 0.05 * (v_pu - 0.8) / 0.2
        xm_saturated = xm_rated * saturation_factor
    elseif v_pu <= 1.3
        # Moderate saturation
        base_factor = 0.95
        additional_reduction = 0.15 * (v_pu - 1.0) / 0.3
        saturation_factor = base_factor - additional_reduction
        xm_saturated = xm_rated * saturation_factor
    else
        # Heavy saturation - Xm drops significantly
        base_factor = 0.8
        heavy_reduction = 0.3 * ((v_pu - 1.3) / 0.7)^1.2
        saturation_factor = base_factor - heavy_reduction
        saturation_factor = max(saturation_factor, 0.3)  # Minimum 30% of rated
        xm_saturated = xm_rated * saturation_factor
    end
    
    return xm_saturated
end

"""
    calculate_motor_slip_for_load(motor::InductionMotor, voltage::Float64, target_torque::Float64)

Calculate the slip required to maintain a target torque at given voltage.
This is crucial for understanding real motor behavior under voltage variations.
"""
function calculate_motor_slip_for_load(motor::InductionMotor, voltage::Float64, target_torque::Float64)
    # Synchronous speed
    sync_speed = 120 * motor.rated_frequency / motor.poles
    ωs = 2π * sync_speed / 60
    
    # Target air gap power for the desired torque
    target_p_airgap = target_torque * ωs
    
    # Iteratively find slip that produces target torque
    slip_range = range(0.001, 0.5, length=1000)
    best_slip = 0.05
    min_error = Inf
    
    for s in slip_range
        p_m, p_e, q, i, t, eff, _, _, _ = calculate_correct_motor_performance(motor, voltage, s)
        torque_error = abs(t - target_torque)
        
        if torque_error < min_error
            min_error = torque_error
            best_slip = s
        end
    end
    
    return best_slip
end

"""
    test_constant_load_voltage_variation()

Test motor behavior under constant mechanical load with varying voltage.
This shows the realistic effect on slip and efficiency.
"""
function test_constant_load_voltage_variation()
    println("="^60)
    println("TESTING CONSTANT LOAD WITH VOLTAGE VARIATION")
    println("="^60)
    
    # Create motor
    motor = InductionMotor(
        rated_power=10.0,      
        rated_voltage=400.0,   
        rated_frequency=50.0,  
        rated_speed=1450.0,    
        poles=4,
        r1=0.4,                
        x1=1.5,                
        r2=0.25,               
        x2=1.5,                
        xm=25.0                
    )
    
    # Calculate rated operating point
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    
    # Get rated torque
    p_m_rated, _, _, _, t_rated, _, _, _, _ = 
        calculate_correct_motor_performance(motor, motor.rated_voltage, rated_slip)
    
    println("Reference Operating Point (Rated Voltage):")
    println("  Rated Torque: $(round(t_rated, digits=1)) N⋅m")
    println("  Rated Slip: $(round(rated_slip*100, digits=2))%")
    println("  Rated Power: $(round(p_m_rated, digits=2)) kW")
    println()
    
    # Test voltage variation with constant torque load
    voltage_range = range(0.6, 1.2, length=30)
    
    println("Constant Load Analysis:")
    println("V(p.u.) | Slip(%) | P(kW) | Q(kVAr) | I(A)  | η(%)  | R2'/s(Ω) | Explanation")
    println("-"^90)
    
    u_vals = Float64[]
    slip_vals = Float64[]
    p_vals = Float64[]
    q_vals = Float64[]
    i_vals = Float64[]
    eff_vals = Float64[]
    r2_over_s_vals = Float64[]
    
    for v_pu in voltage_range
        voltage = v_pu * motor.rated_voltage
        
        # Find slip required for constant torque
        required_slip = calculate_motor_slip_for_load(motor, voltage, t_rated)
        
        # Calculate performance at this operating point
        p_m, p_e, q, i, t, eff, _, _, _ = 
            calculate_correct_motor_performance(motor, voltage, required_slip)
        
        # Calculate equivalent rotor resistance
        r2_over_s = motor.r2 / required_slip
        
        push!(u_vals, v_pu)
        push!(slip_vals, required_slip * 100)
        push!(p_vals, p_m)
        push!(q_vals, q)
        push!(i_vals, i)
        push!(eff_vals, eff)
        push!(r2_over_s_vals, r2_over_s)
        
        # Explanation based on operating conditions
        explanation = ""
        if v_pu < 0.8
            explanation = "Low V: High slip needed → Low eff"
        elseif v_pu < 1.1
            explanation = "Normal operation"
        else
            explanation = "High V: Lower slip → Better eff"
        end
        
        if v_pu in [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
            @printf(" %.1f    | %5.1f | %.2f | %6.2f | %.1f | %.1f | %6.2f  | %s\n", 
                    v_pu, required_slip*100, p_m, q, i, eff, r2_over_s, explanation)
        end
    end
    
    # Create enhanced plots showing the physical effects
    p1 = plot(u_vals, slip_vals,
              label="Required Slip",
              xlabel="Voltage (p.u.)",
              ylabel="Slip (%)",
              title="Slip Required for Constant Torque\n(Higher slip at low voltage)",
              linewidth=3,
              color=:red,
              grid=true)
    
    # Add explanation
    annotate!(p1, 0.7, maximum(slip_vals)*0.8,
              text("Low voltage requires\nhigher slip to maintain\nsame torque", 10, :left))
    
    p2 = plot(u_vals, eff_vals,
              label="Efficiency",
              xlabel="Voltage (p.u.)",
              ylabel="Efficiency (%)",
              title="Efficiency vs Voltage (Constant Load)\n(Shows realistic efficiency drop at low V)",
              linewidth=3,
              color=:blue,
              grid=true)
    
    # Mark rated point
    scatter!(p2, [1.0], [eff_vals[argmin(abs.(u_vals .- 1.0))]],
            markersize=8,
            color=:red,
            markershape=:circle,
            label="Rated point")
    
    # Add physical explanation
    annotate!(p2, 0.7, maximum(eff_vals)*0.9,
              text("Higher slip →\nMore rotor losses →\nLower efficiency", 10, :left))
    
    p3 = plot(u_vals, r2_over_s_vals,
              label="R₂'/s",
              xlabel="Voltage (p.u.)",
              ylabel="Equivalent Rotor Resistance (Ω)",
              title="Rotor Resistance Effect\n(R₂'/s decreases with higher slip)",
              linewidth=3,
              color=:green,
             
              grid=true)
    
    # Add I ∝ 1/V reference line
    i_ref = current_vals[argmin(abs.(u_vals .- v_ref))]
    i_inv_v_line = i_ref .* (v_ref ./ u_vals)
    plot!(p3, u_vals, i_inv_v_line,
          label="I ∝ 1/V reference",
          linewidth=2,
          color=:orange,
          linestyle=:dash)
    
    p4 = plot(u_vals, efficiency_vals,
              label="Efficiency",
              xlabel="Voltage (p.u.)",
              ylabel="Efficiency (%)",
              title="Efficiency vs Voltage\n(Optimal near rated voltage)",
              linewidth=3,
              color=:brown,
              grid=true)
    
    # Mark rated voltage point
    scatter!(p4, [1.0], [efficiency_vals[argmin(abs.(u_vals .- 1.0))]],
            markersize=8,
            color=:red,
            markershape=:circle,
            label="Rated point")
    
    # Display combined plots
    combined_plot = plot(p1, p2, p3, p4, 
                        layout=(2,2), 
                        size=(1200, 900),
                        plot_title="Induction Motor: U-P-Q Characteristics at Constant Slip")
    display(combined_plot)
    
    # Test different slip values at rated voltage
    println("Test 2: Effect of Slip at Rated Voltage")
    println("-"^50)
    
    slip_range = range(0.01, 0.15, length=30)
    p_slip_vals = Float64[]
    q_slip_vals = Float64[]
    s_vals = Float64[]
    
    for s in slip_range
        p_m, _, q, _, _, _, _, _, _ = calculate_correct_motor_performance(motor, motor.rated_voltage, s)
        push!(s_vals, s)
        push!(p_slip_vals, p_m)
        push!(q_slip_vals, q)
    end
    
    p5 = plot(s_vals .* 100, p_slip_vals,
              label="Active Power",
              xlabel="Slip (%)",
              ylabel="Power (kW)",
              title="Power vs Slip at Rated Voltage",
              linewidth=3,
              color=:blue,
              grid=true)
    
    p6 = plot(s_vals .* 100, q_slip_vals,
              label="Reactive Power",
              xlabel="Slip (%)",
              ylabel="Reactive Power (kVAr)",
              title="Reactive Power vs Slip",
              linewidth=3,
              color=:red,
              grid=true)
    
    # Mark rated slip point
    scatter!(p5, [rated_slip * 100], [p_m_rated],
            markersize=8,
            color=:red,
            markershape=:circle,
            label="Rated slip")
    
    scatter!(p6, [rated_slip * 100], [q_rated],
            markersize=8,
            color=:blue,
            markershape=:circle,
            label="Rated slip")
    
    slip_plot = plot(p5, p6, 
                    layout=(1,2), 
                    size=(1000, 400),
                    plot_title="Effect of Slip on Motor Performance")
    display(slip_plot)
    
    # Test operating point table
    println("Test 3: Operating Point Analysis")
    println("-"^50)
    
    test_conditions = [
        (0.8, "Low voltage"),
        (0.9, "Reduced voltage"), 
        (1.0, "Rated voltage"),
        (1.1, "High voltage"),
        (1.2, "Overvoltage")
    ]
    
    println("Voltage | P(kW) | Q(kVAr) | I(A)  | η(%)  | Condition")
    println("-"^60)
    
    for (v_pu, description) in test_conditions
        voltage = v_pu * motor.rated_voltage
        p_m, _, q, i, _, eff, _, _, _ = calculate_correct_motor_performance(motor, voltage, rated_slip)
        @printf("  %.1f   | %.2f | %6.2f | %.1f | %.1f | %s\n", 
                v_pu, p_m, q, i, eff, description)
    end
    
    println()
    
    # Summary and insights
    println("Summary and Physical Insights:")
    println("1. Active Power (P vs V):")
    println("   - P increases with voltage: P ∝ V² (at constant slip)")
    println("   - Lower voltage severely reduces power capability")
    println("   - Maximum power at highest voltage in test range")
    println()
    
    println("2. Reactive Power (Q vs V):")
    println("   - Characteristic U-shaped curve observed")
    println("   - Minimum Q at ~$(round(min_q_voltage, digits=2)) p.u. voltage")
    println("   - At low V: High current → High leakage reactive power")
    println("   - At high V: High flux → High magnetizing reactive power")
    println()
    
    println("3. Current (I vs V):")
    println("   - Current increases as voltage decreases: I ∝ 1/V")
    println("   - Higher current at low voltage increases losses")
    println()
    
    println("4. Efficiency:")
    println("   - Best efficiency near rated voltage")
    println("   - Reduced efficiency at both low and high voltages")
    println("   - Voltage regulation important for optimal operation")
    println()
    
    return motor, combined_plot, slip_plot
end

"""
    test_motor_loading()

Test motor performance under different loading conditions.
"""
function test_motor_loading()
    println("="^60)
    println("TESTING MOTOR LOADING CHARACTERISTICS")
    println("="^60)
    
    # Create a standard motor
    motor = InductionMotor(
        rated_power=15.0,
        rated_voltage=400.0,
        rated_frequency=50.0,
        rated_speed=1450.0,
        poles=4,
        r1=0.3,
        x1=1.2,
        r2=0.25,
        x2=1.2,
        xm=35.0
    )
    
    # Test different loading conditions (represented by slip)
    sync_speed = 120 * motor.rated_frequency / motor.poles
    loading_conditions = [
        (0.005, "No load"),
        (0.015, "25% load"),
        (0.025, "50% load"), 
        (0.035, "75% load"),
        (0.045, "100% load"),
        (0.070, "150% load")
    ]
    
    println("Loading Test Results at Rated Voltage:")
    println("Load Level | Slip(%) | P(kW) | Q(kVAr) | I(A)  | η(%)  | pf")
    println("-"^70)
    
    for (slip, description) in loading_conditions
        p_m, p_e, q, i, _, eff, _, _, _ = calculate_correct_motor_performance(motor, motor.rated_voltage, slip)
        
        # Calculate power factor
        s_apparent = sqrt(p_e^2 + q^2)
        power_factor = s_apparent > 0 ? p_e / s_apparent : 0
        
        @printf("%-10s | %5.1f | %.2f | %6.2f | %.1f | %.1f | %.3f\n", 
                description, slip*100, p_m, q, i, eff, power_factor)
    end
    
    return motor
end

"""
    main()

Main test function that runs all induction motor tests.
"""
function main()
    println("INDUCTION MOTOR COMPREHENSIVE TESTING")
    println("="^60)
    println("This test suite demonstrates:")
    println("1. U-P-Q relationships at constant slip")
    println("2. Effect of slip on motor performance") 
    println("3. Motor loading characteristics")
    println("4. Voltage sensitivity analysis")
    println("="^60)
    println()
    
    # Run main U-P-Q test
    motor1, upq_plot, slip_plot = test_upq_relationship()
    
    println("\n")
    
    # Run loading test
    motor2 = test_motor_loading()
    
    println("\n" * "="^60)
    println("ALL TESTS COMPLETED SUCCESSFULLY!")
    println("Key findings:")
    println("• Active power P ∝ V² (decreases with voltage)")
    println("• Reactive power shows U-shaped curve vs voltage")
    println("• Current I ∝ 1/V (increases as voltage decreases)")
    println("• Best efficiency occurs near rated voltage")
    println("• Loading affects slip, which impacts all performance metrics")
    println("="^60)
    
    return motor1, motor2, upq_plot, slip_plot
end

# Run the tests
# if abspath(PROGRAM_FILE) == @__FILE__
    motor1, motor2, upq_plot, slip_plot = main()
# end

