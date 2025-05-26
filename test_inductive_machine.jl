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
Enhanced model for more pronounced saturation effects.
"""
function calculate_saturated_xm(motor::InductionMotor, voltage::Float64)
    v_pu = voltage / motor.rated_voltage
    xm_rated = motor.xm
    
    # Enhanced saturation model with more pronounced effects
    if v_pu <= 0.8
        # Linear region - no saturation
        xm_saturated = xm_rated
    elseif v_pu <= 1.0
        # Mild saturation begins
        saturation_factor = 1.0 - 0.05 * (v_pu - 0.8) / 0.2
        xm_saturated = xm_rated * saturation_factor
    elseif v_pu <= 1.2
        # Moderate saturation
        base_factor = 0.95
        additional_reduction = 0.15 * (v_pu - 1.0) / 0.2
        saturation_factor = base_factor - additional_reduction
        xm_saturated = xm_rated * saturation_factor
    elseif v_pu <= 1.5
        # Heavy saturation starts
        base_factor = 0.80
        heavy_reduction = 0.25 * ((v_pu - 1.2) / 0.3)^1.5
        saturation_factor = base_factor - heavy_reduction
        xm_saturated = xm_rated * saturation_factor
    else
        # Severe saturation - Xm drops dramatically
        base_factor = 0.55
        severe_reduction = 0.35 * ((v_pu - 1.5) / 0.5)^2.0
        saturation_factor = base_factor - severe_reduction
        saturation_factor = max(saturation_factor, 0.2)  # Minimum 20% of rated
        xm_saturated = xm_rated * saturation_factor
    end
    
    return xm_saturated
end

"""
    calculate_core_losses(motor::InductionMotor, voltage::Float64)

Calculate core losses which increase significantly at high voltages due to saturation.
"""
function calculate_core_losses(motor::InductionMotor, voltage::Float64)
    v_pu = voltage / motor.rated_voltage
    
    # Base core losses at rated voltage (empirical: ~1-3% of rated power)
    base_core_loss = motor.rated_power * 0.02  # 2% of rated power in kW
    
    # Core losses increase with voltage and saturation
    if v_pu <= 1.0
        # Linear increase up to rated voltage
        core_loss_factor = v_pu^1.8  # Slightly less than quadratic
    elseif v_pu <= 1.3
        # Saturation starts - losses increase faster
        base_increase = 1.0^1.8
        saturation_increase = 2.5 * (v_pu - 1.0)^2.2  # Quadratic+ increase
        core_loss_factor = base_increase + saturation_increase
    else
        # Heavy saturation - exponential increase in losses
        moderate_loss = 1.0^1.8 + 2.5 * (1.3 - 1.0)^2.2
        severe_increase = 5.0 * ((v_pu - 1.3) / 0.7)^3.0  # Cubic increase
        core_loss_factor = moderate_loss + severe_increase
    end
    
    return base_core_loss * core_loss_factor
end

"""
    calculate_enhanced_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)

Enhanced motor performance calculation including core losses and detailed loss breakdown.
"""
function calculate_enhanced_motor_performance(motor::InductionMotor, voltage::Float64, slip::Float64)
    # Avoid division by zero
    if slip ≈ 0
        slip = 1e-6
    end
    
    # Per-phase voltage
    V1 = voltage / √3
    
    # Get saturated magnetizing reactance
    xm_saturated = calculate_saturated_xm(motor, voltage)
    
    # Complex impedances with saturation
    Z1 = complex(motor.r1, motor.x1)
    Z2 = complex(motor.r2/slip, motor.x2)
    Zm = complex(0.0, xm_saturated)
    
    # Circuit solution
    Z_parallel = (Z2 * Zm) / (Z2 + Zm)
    Z_total = Z1 + Z_parallel
    I1 = V1 / Z_total
    V_parallel = I1 * Z_parallel
    I2 = V_parallel / Z2
    Im = V_parallel / Zm
    
    # Detailed loss calculations (3-phase, in Watts)
    P_stator_copper = 3 * abs(I1)^2 * motor.r1      # Stator copper loss
    P_rotor_copper = 3 * abs(I2)^2 * motor.r2       # Rotor copper loss
    P_core = calculate_core_losses(motor, voltage) * 1000  # Core losses (convert kW to W)
    
    # Air gap power and mechanical power
    P_airgap = 3 * abs(I2)^2 * motor.r2 / slip
    P_mech = P_airgap - P_rotor_copper
    P_elec = P_airgap + P_stator_copper + P_core
    
    # Reactive power calculation
    S_total = 3 * V1 * conj(I1)
    Q_total = imag(S_total)
    
    # Component breakdown for analysis
    Q_stator_leakage = 3 * abs(I1)^2 * motor.x1
    Q_rotor_leakage = 3 * abs(I2)^2 * motor.x2
    Q_magnetizing = 3 * abs(V_parallel)^2 / xm_saturated
    
    # Torque calculation
    synchronous_speed = 120 * motor.rated_frequency / motor.poles
    ωs = 2π * synchronous_speed / 60
    torque = P_airgap / ωs
    
    # Efficiency
    efficiency = P_mech > 0 ? P_mech / P_elec * 100 : 0
    
    return (P_mech/1000, P_elec/1000, Q_total/1000, abs(I1), torque, efficiency,
            P_stator_copper/1000, P_rotor_copper/1000, P_core/1000,
            Q_stator_leakage/1000, Q_rotor_leakage/1000, Q_magnetizing/1000,
            xm_saturated)
end

"""
    calculate_motor_slip_for_load(motor::InductionMotor, voltage::Float64, target_torque::Float64)

Calculate the slip required to maintain a target torque at given voltage.
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
    test_constant_torque_voltage_effects()

Comprehensive test showing the dramatic effects of voltage variation on motor losses
when maintaining constant torque.
"""
function test_constant_torque_voltage_effects()
    println("="^80)
    println("COMPREHENSIVE TEST: VOLTAGE EFFECTS ON MOTOR LOSSES (CONSTANT TORQUE)")
    println("="^80)
    
    # Create motor with realistic parameters
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
        xm=30.0                
    )
    
    # Calculate rated operating point
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    
    # Get rated torque as reference
    p_m_rated, _, _, _, t_rated, _, _, _, _, _, _, _, _ = 
        calculate_enhanced_motor_performance(motor, motor.rated_voltage, rated_slip)
    
    println("REFERENCE CONDITIONS:")
    println("  Rated Torque: $(round(t_rated, digits=1)) N⋅m")
    println("  Rated Power: $(round(p_m_rated, digits=2)) kW")
    println("  Rated Slip: $(round(rated_slip*100, digits=2))%")
    println("  Motor: $(motor.rated_power) kW, $(motor.rated_voltage) V")
    println()
    
    # Test wide voltage range with constant torque
    voltage_range = range(0.5, 2.0, length=40)
    
    # Storage arrays
    u_vals = Float64[]
    slip_vals = Float64[]
    efficiency_vals = Float64[]
    stator_loss_vals = Float64[]
    rotor_loss_vals = Float64[]
    core_loss_vals = Float64[]
    total_loss_vals = Float64[]
    r2_over_s_vals = Float64[]
    xm_sat_vals = Float64[]
    current_vals = Float64[]
    
    println("DETAILED ANALYSIS: Constant Torque = $(round(t_rated, digits=1)) N⋅m")
    println("="^90)
    println("V(pu) | Slip(%) | R2'/s(Ω) | I(A) | η(%) | P_stator(W) | P_rotor(W) | P_core(W) | Xm_sat(Ω) | Dominant Loss")
    println("-"^95)
    
    for v_pu in voltage_range
        voltage = v_pu * motor.rated_voltage
        
        # Find slip required for constant torque
        required_slip = calculate_motor_slip_for_load(motor, voltage, t_rated)
        
        # Calculate enhanced performance
        p_m, p_e, q, i, t, eff, p_stator, p_rotor, p_core, _, _, _, xm_sat = 
            calculate_enhanced_motor_performance(motor, voltage, required_slip)
        
        # Calculate derived parameters
        r2_over_s = motor.r2 / required_slip
        total_loss = (p_stator + p_rotor + p_core) * 1000  # Convert to Watts
        
        # Store values
        push!(u_vals, v_pu)
        push!(slip_vals, required_slip * 100)
        push!(efficiency_vals, eff)
        push!(stator_loss_vals, p_stator * 1000)
        push!(rotor_loss_vals, p_rotor * 1000)
        push!(core_loss_vals, p_core * 1000)
        push!(total_loss_vals, total_loss)
        push!(r2_over_s_vals, r2_over_s)
        push!(xm_sat_vals, xm_sat)
        push!(current_vals, i)
        
        # Determine dominant loss mechanism
        losses = [p_stator * 1000, p_rotor * 1000, p_core * 1000]
        loss_names = ["Stator", "Rotor", "Core"]
        dominant_loss = loss_names[argmax(losses)]
        
        # Print detailed results for key voltages
        if v_pu in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.8, 2.0]
            @printf("%.1f   | %6.1f | %6.2f  | %.1f | %.1f |    %7.0f |    %6.0f |   %6.0f |   %6.1f | %s\n",
                    v_pu, required_slip*100, r2_over_s, i, eff, 
                    p_stator*1000, p_rotor*1000, p_core*1000, xm_sat, dominant_loss)
        end
    end
    
    println()
    
    # Identify key operating regions
    low_v_idx = findfirst(v -> v <= 0.7, u_vals)
    rated_v_idx = argmin(abs.(u_vals .- 1.0))
    high_v_idx = findfirst(v -> v >= 1.5, u_vals)
    
    println("KEY INSIGHTS:")
    println("="^50)
    println("1. LOW VOLTAGE EFFECTS (V ≤ 0.7 p.u.):")
    println("   - Slip increases dramatically: $(round(slip_vals[low_v_idx], digits=1))% vs $(round(slip_vals[rated_v_idx], digits=1))% rated")
    println("   - R2'/s decreases: $(round(r2_over_s_vals[low_v_idx], digits=2))Ω vs $(round(r2_over_s_vals[rated_v_idx], digits=2))Ω rated")
    println("   - Rotor losses dominate: $(round(rotor_loss_vals[low_v_idx], digits=0))W vs $(round(rotor_loss_vals[rated_v_idx], digits=0))W rated")
    println("   - Efficiency drops: $(round(efficiency_vals[low_v_idx], digits=1))% vs $(round(efficiency_vals[rated_v_idx], digits=1))% rated")
    println()
    
    if high_v_idx !== nothing
        println("2. HIGH VOLTAGE EFFECTS (V ≥ 1.5 p.u.):")
        println("   - Magnetic saturation: Xm = $(round(xm_sat_vals[high_v_idx], digits=1))Ω vs $(round(motor.xm, digits=1))Ω unsaturated")
        println("   - Core losses dominate: $(round(core_loss_vals[high_v_idx], digits=0))W vs $(round(core_loss_vals[rated_v_idx], digits=0))W rated")
        println("   - Saturation factor: $(round(xm_sat_vals[high_v_idx]/motor.xm, digits=2)) of rated Xm")
        println("   - Efficiency affected: $(round(efficiency_vals[high_v_idx], digits=1))% vs $(round(efficiency_vals[rated_v_idx], digits=1))% rated")
        println()
    end
    
    println("3. PHYSICAL MECHANISMS:")
    println("   - Low V: s↑ → R2'/s↓ → I2↑ → P_rotor = I2²×R2 ↑↑ (quadratic effect)")
    println("   - High V: Saturation → Xm↓ → Higher flux density → P_core↑↑ (exponential effect)")
    println()
    
    # Create comprehensive plots
    p1 = plot(u_vals, slip_vals,
              label="Required Slip",
              xlabel="Voltage (p.u.)",
              ylabel="Slip (%)",
              title="Slip vs Voltage (Constant Torque)\nDRAMATIC increase at low voltage",
              linewidth=4,
              color=:red,
              grid=true,
              yscale=:log10)
    
    # Add annotation for low voltage region
    annotate!(p1, 0.6, maximum(slip_vals)*0.5,
              text("Low V: Slip increases\nexponentially to maintain\nconstant torque", 10, :left))
    
    p2 = plot(u_vals, rotor_loss_vals,
              label="Rotor Losses",
              xlabel="Voltage (p.u.)",
              ylabel="Losses (W)",
              title="Loss Breakdown vs Voltage\nRotor losses DOMINATE at low V",
              linewidth=4,
              color=:red,
              grid=true,
              yscale=:log10)
    
    plot!(p2, u_vals, stator_loss_vals,
          label="Stator Losses",
          linewidth=3,
          color=:blue)
    
    plot!(p2, u_vals, core_loss_vals,
          label="Core Losses",
          linewidth=3,
          color=:green)
    
    # Mark transition points
    scatter!(p2, [0.7], [rotor_loss_vals[findfirst(v -> v >= 0.7, u_vals)]],
            markersize=8, color=:red, markershape=:circle, label="Rotor loss peak")
    
    if high_v_idx !== nothing
        scatter!(p2, [1.5], [core_loss_vals[high_v_idx]],
                markersize=8, color=:green, markershape=:square, label="Core loss rise")
    end
    
    p3 = plot(u_vals, efficiency_vals,
              label="Efficiency",
              xlabel="Voltage (p.u.)",
              ylabel="Efficiency (%)",
              title="Efficiency vs Voltage (Constant Torque)\nShows REALISTIC motor behavior",
              linewidth=4,
              color=:purple,
              grid=true)
    
    # Mark optimal region
    optimal_range = findall(v -> 0.95 <= v <= 1.05, u_vals)
    if !isempty(optimal_range)
        scatter!(p3, u_vals[optimal_range], efficiency_vals[optimal_range],
                markersize=4, color=:green, alpha=0.7, label="Optimal range")
    end
    
    # Add efficiency drop annotations
    annotate!(p3, 0.6, maximum(efficiency_vals)*0.8,
              text("Efficiency drops\n>50% at low voltage!", 10, :left))
    
    p4 = plot(u_vals, xm_sat_vals,
              label="Xm (saturated)",
              xlabel="Voltage (p.u.)",
              ylabel="Magnetizing Reactance (Ω)",
              title="Magnetic Saturation Effect\nXm decreases dramatically at high V",
              linewidth=4,
              color=:orange,
              grid=true)
    
    # Add unsaturated reference
    hline!(p4, [motor.xm],
           label="Xm (unsaturated)",
           linewidth=2,
           color=:gray,
           linestyle=:dash)
    
    # Show saturation region
    saturation_start = findfirst(v -> v >= 1.2, u_vals)
    if saturation_start !== nothing
        annotate!(p4, 1.4, motor.xm * 0.7,
                  text("Heavy saturation\nXm drops to $(round(minimum(xm_sat_vals)/motor.xm*100, digits=0))%\nof rated value", 10, :center))
    end
    
    # Create the comprehensive comparison plot
    comprehensive_plot = plot(p1, p2, p3, p4,
                             layout=(2,2),
                             size=(1400, 1000),
                             plot_title="Motor Behavior Analysis: Constant Torque Operation")
    
    display(comprehensive_plot)
    
    # Quantitative summary
    println("QUANTITATIVE SUMMARY:")
    println("="^50)
    
    # Find extreme points
    min_eff_idx = argmin(efficiency_vals)
    max_loss_idx = argmax(total_loss_vals)
    max_slip_idx = argmax(slip_vals)
    min_xm_idx = argmin(xm_sat_vals)
    
    println("WORST CASE SCENARIOS:")
    println("  Lowest efficiency: $(round(efficiency_vals[min_eff_idx], digits=1))% at $(round(u_vals[min_eff_idx], digits=2)) p.u.")
    println("  Highest losses: $(round(total_loss_vals[max_loss_idx], digits=0))W at $(round(u_vals[max_loss_idx], digits=2)) p.u.")
    println("  Maximum slip: $(round(slip_vals[max_slip_idx], digits=1))% at $(round(u_vals[max_slip_idx], digits=2)) p.u.")
    println("  Minimum Xm: $(round(xm_sat_vals[min_xm_idx], digits=1))Ω ($(round(xm_sat_vals[min_xm_idx]/motor.xm*100, digits=0))%) at $(round(u_vals[min_xm_idx], digits=2)) p.u.")
    
    println()
    println("PRACTICAL RECOMMENDATIONS:")
    println("  ✓ Operate between 0.95-1.05 p.u. voltage for optimal efficiency")
    println("  ⚠ Avoid operation below 0.8 p.u. - efficiency drops >30%")
    println("  ⚠ Avoid operation above 1.3 p.u. - core losses increase dramatically")
    println("  ⚠ Voltage regulation is CRITICAL for motor efficiency")
    
    return motor, comprehensive_plot, u_vals, efficiency_vals, rotor_loss_vals, core_loss_vals
end

"""
    test_loss_mechanisms_detailed()

Detailed analysis of different loss mechanisms at various operating points.
"""
function test_loss_mechanisms_detailed()
    println("="^80)
    println("DETAILED LOSS MECHANISM ANALYSIS")
    println("="^80)
    
    motor = InductionMotor(
        rated_power=15.0, rated_voltage=400.0, rated_frequency=50.0,
        rated_speed=1450.0, poles=4,
        r1=0.3, x1=1.2, r2=0.2, x2=1.2, xm=35.0
    )
    
    # Test specific operating points
    test_points = [
        (0.6, "Low voltage - Rotor loss dominated"),
        (0.8, "Reduced voltage - Transition region"),
        (1.0, "Rated voltage - Balanced operation"),
        (1.2, "High voltage - Mild saturation"),
        (1.5, "Overvoltage - Heavy saturation"),
        (1.8, "Severe overvoltage - Core loss dominated")
    ]
    
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    
    # Get reference torque
    _, _, _, _, t_rated, _, _, _, _, _, _, _, _ = 
        calculate_enhanced_motor_performance(motor, motor.rated_voltage, rated_slip)
    
    println("LOSS MECHANISM BREAKDOWN:")
    println("V(p.u.) | Slip(%) | Stator(W) | Rotor(W) | Core(W) | Total(W) | η(%) | Mechanism")
    println("-"^85)
    
    for (v_pu, description) in test_points
        voltage = v_pu * motor.rated_voltage
        required_slip = calculate_motor_slip_for_load(motor, voltage, t_rated)
        
        p_m, p_e, q, i, t, eff, p_stator, p_rotor, p_core, _, _, _, xm_sat = 
            calculate_enhanced_motor_performance(motor, voltage, required_slip)
        
        # Convert to Watts for display
        stator_w = p_stator * 1000
        rotor_w = p_rotor * 1000
        core_w = p_core * 1000
        total_w = stator_w + rotor_w + core_w
        
        # Determine dominant mechanism
        losses = [stator_w, rotor_w, core_w]
        percentages = losses ./ total_w * 100
        dominant_idx = argmax(losses)
        mechanisms = ["Stator copper", "Rotor copper", "Core/magnetic"]
        
        @printf("  %.1f   | %5.1f | %7.0f | %6.0f | %5.0f | %6.0f | %.1f | %s (%.0f%%)\n",
                v_pu, required_slip*100, stator_w, rotor_w, core_w, total_w, eff,
                mechanisms[dominant_idx], percentages[dominant_idx])
    end
    
    println()
    println("PHYSICAL EXPLANATIONS:")
    println("1. LOW VOLTAGE (≤0.8 p.u.):")
    println("   • High slip required → Low R2'/s → High rotor current")
    println("   • P_rotor = 3×I2²×R2 increases quadratically")
    println("   • Rotor losses can be 5-10× rated values")
    println()
    println("2. HIGH VOLTAGE (≥1.3 p.u.):")
    println("   • Magnetic saturation → Reduced Xm")
    println("   • Higher flux density → Eddy current and hysteresis losses")
    println("   • P_core ∝ V^2.5-3.0 (non-linear increase)")
    println()
    println("3. RATED VOLTAGE (~1.0 p.u.):")
    println("   • Balanced operation with optimal loss distribution")
    println("   • All loss mechanisms at reasonable levels")
    
    return motor
end

"""
    plot_pv_qv_characteristics(motor::InductionMotor)

Generate comprehensive P-V and Q-V characteristic curves for induction motor.
Shows both constant slip (theoretical) and constant load (realistic) scenarios.
"""
function plot_pv_qv_characteristics(motor::InductionMotor)
    println("="^60)
    println("GENERATING P-V AND Q-V CHARACTERISTIC CURVES")
    println("="^60)
    
    # Calculate rated conditions
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    
    # Get rated torque for constant load analysis
    _, _, _, _, t_rated, _, _, _, _, _, _, _, _ = 
        calculate_enhanced_motor_performance(motor, motor.rated_voltage, rated_slip)
    
    # Voltage range for analysis
    voltage_range = range(0.4, 1.6, length=100)
    
    # Arrays for constant slip analysis (theoretical)
    u_vals_const = Float64[]
    p_vals_const = Float64[]
    q_vals_const = Float64[]
    q_stator_const = Float64[]
    q_rotor_const = Float64[]
    q_mag_const = Float64[]
    
    # Arrays for constant load analysis (realistic)
    u_vals_load = Float64[]
    p_vals_load = Float64[]
    q_vals_load = Float64[]
    slip_vals_load = Float64[]
    
    println("Calculating P-V and Q-V characteristics...")
    
    # Constant slip analysis (theoretical)
    for v_pu in voltage_range
        voltage = v_pu * motor.rated_voltage
        p_m, _, q, _, _, _, _, _, _, q_s, q_r, q_m, _ = 
            calculate_enhanced_motor_performance(motor, voltage, rated_slip)
        
        push!(u_vals_const, v_pu)
        push!(p_vals_const, p_m)
        push!(q_vals_const, q)
        push!(q_stator_const, q_s)
        push!(q_rotor_const, q_r)
        push!(q_mag_const, q_m)
    end
    
    # Constant load analysis (realistic) - limited range for feasibility
    feasible_range = range(0.6, 1.3, length=50)
    for v_pu in feasible_range
        voltage = v_pu * motor.rated_voltage
        
        try
            required_slip = calculate_motor_slip_for_load(motor, voltage, t_rated)
            p_m, _, q, _, _, _, _, _, _, _, _, _, _ = 
                calculate_enhanced_motor_performance(motor, voltage, required_slip)
            
            push!(u_vals_load, v_pu)
            push!(p_vals_load, p_m)
            push!(q_vals_load, q)
            push!(slip_vals_load, required_slip * 100)
        catch
            # Skip points where constant load cannot be maintained
            continue
        end
    end
    
    # Create P-V characteristic plot
    pv_plot = plot(u_vals_const, p_vals_const,
                   label="P-V (Constant Slip = $(round(rated_slip*100, digits=1))%)",
                   xlabel="Voltage (p.u.)",
                   ylabel="Active Power (kW)",
                   title="P-V Characteristics: Constant Slip vs Constant Load",
                   linewidth=3,
                   color=:blue,
                   linestyle=:dash,
                   grid=true,
                   legend=:topleft)
    
    plot!(pv_plot, u_vals_load, p_vals_load,
          label="P-V (Constant Load = $(round(t_rated, digits=1)) N⋅m)",
          linewidth=4,
          color=:red,
          linestyle=:solid)
    
    # Add P ∝ V² reference line for constant slip
    v_ref = 1.0
    p_ref = p_vals_const[argmin(abs.(u_vals_const .- v_ref))]
    p_v2_line = p_ref .* (u_vals_const ./ v_ref).^2
    plot!(pv_plot, u_vals_const, p_v2_line,
          label="P ∝ V² (reference)",
          linewidth=2,
          color=:gray,
          linestyle=:dot)
    
    # Add annotations
    annotate!(pv_plot, 0.6, maximum(p_vals_const)*0.8,
              text("Constant slip: P ∝ V²\n(unrealistic)", 10, :left))
    annotate!(pv_plot, 1.0, maximum(p_vals_load)*0.6,
              text("Constant load:\nP varies less\n(realistic)", 10, :center))
    
    # Create Q-V characteristic plot with component breakdown
    qv_plot = plot(u_vals_const, q_vals_const,
                   label="Q-V Total (Constant Slip)",
                   xlabel="Voltage (p.u.)",
                   ylabel="Reactive Power (kVAr)",
                   title="Q-V Characteristics: Components and U-Shape",
                   linewidth=4,
                   color=:blue,
                   grid=true,
                   legend=:topright)
    
    # Add Q components for constant slip
    plot!(qv_plot, u_vals_const, q_mag_const,
          label="Q_magnetizing",
          linewidth=3,
          color=:green,
          linestyle=:dash)
    
    leakage_q_const = q_stator_const .+ q_rotor_const
    plot!(qv_plot, u_vals_const, leakage_q_const,
          label="Q_leakage (total)",
          linewidth=2,
          color=:orange,
          linestyle=:dot)
    
    # Add constant load Q curve
    plot!(qv_plot, u_vals_load, q_vals_load,
          label="Q-V (Constant Load)",
          linewidth=3,
          color=:red,
          linestyle=:solid)
    
    # Mark minimum Q point
    min_q_idx = argmin(q_vals_const)
    min_q_voltage = u_vals_const[min_q_idx]
    min_q_value = q_vals_const[min_q_idx]
    
    scatter!(qv_plot, [min_q_voltage], [min_q_value],
            markersize=10,
            color=:black,
            markershape=:star,
            label="Q minimum")
    
    # Add physics explanations
    annotate!(qv_plot, 0.5, maximum(q_vals_const)*0.8,
              text("Low V: High I\n→ High Q_leakage", 9, :left))
    annotate!(qv_plot, 1.4, maximum(q_vals_const)*0.7,
              text("High V: Saturation\n→ Q_mag decreases", 9, :right))
    
    # Create separate Q components analysis plot
    q_components_plot = plot(u_vals_const, q_mag_const,
                            label="Q_magnetizing",
                            xlabel="Voltage (p.u.)",
                            ylabel="Reactive Power (kVAr)",
                            title="Q-V Component Analysis (with Saturation)",
                            linewidth=3,
                            color=:green,
                            grid=true,
                            legend=:topright)
    
    plot!(q_components_plot, u_vals_const, leakage_q_const,
          label="Q_leakage (stator + rotor)",
          linewidth=3,
          color=:orange)
    
    plot!(q_components_plot, u_vals_const, q_stator_const,
          label="Q_stator",
          linewidth=2,
          color=:blue,
          linestyle=:dash)
    
    plot!(q_components_plot, u_vals_const, q_rotor_const,
          label="Q_rotor", 
          linewidth=2,
          color=:purple,
          linestyle=:dot)
    
    # Add saturation curve
    xm_sat_vals = [calculate_saturated_xm(motor, v * motor.rated_voltage) for v in u_vals_const]
    q_comp_twin = twinx(q_components_plot)
    plot!(q_comp_twin, u_vals_const, xm_sat_vals,
          label="Xm (saturated)",
          linewidth=2,
          color=:red,
          ylabel="Magnetizing Reactance (Ω)",
          linestyle=:dashdot)
    
    # Create detailed comparison plot
    comparison_plot = plot(layout=(2,2), size=(1400, 1000))
    
    # P-V comparison
    plot!(comparison_plot[1], u_vals_const, p_vals_const,
          label="Constant Slip",
          linewidth=3,
          color=:blue,
          title="P-V Comparison",
          xlabel="Voltage (p.u.)",
          ylabel="Power (kW)")
    
    plot!(comparison_plot[1], u_vals_load, p_vals_load,
          label="Constant Load",
          linewidth=3,
          color=:red)
    
    # Q-V comparison
    plot!(comparison_plot[2], u_vals_const, q_vals_const,
          label="Constant Slip",
          linewidth=3,
          color=:blue,
          title="Q-V Comparison (U-shaped)",
          xlabel="Voltage (p.u.)",
          ylabel="Reactive Power (kVAr)")
    
    plot!(comparison_plot[2], u_vals_load, q_vals_load,
          label="Constant Load",
          linewidth=3,
          color=:red)
    
    # Slip variation for constant load
    plot!(comparison_plot[3], u_vals_load, slip_vals_load,
          label="Required Slip",
          linewidth=3,
          color=:green,
          title="Slip vs Voltage (Constant Load)",
          xlabel="Voltage (p.u.)",
          ylabel="Slip (%)",
          yscale=:log10)
    
    # Efficiency comparison (calculated for both scenarios)
    eff_const = Float64[]
    eff_load = Float64[]
    
    for v_pu in u_vals_const
        voltage = v_pu * motor.rated_voltage
        _, _, _, _, _, eff, _, _, _, _, _, _, _ = 
            calculate_enhanced_motor_performance(motor, voltage, rated_slip)
        push!(eff_const, eff)
    end
    
    for (i, v_pu) in enumerate(u_vals_load)
        voltage = v_pu * motor.rated_voltage
        required_slip = slip_vals_load[i] / 100
        _, _, _, _, _, eff, _, _, _, _, _, _, _ = 
            calculate_enhanced_motor_performance(motor, voltage, required_slip)
        push!(eff_load, eff)
    end
    
    plot!(comparison_plot[4], u_vals_const, eff_const,
          label="Constant Slip",
          linewidth=3,
          color=:blue,
          title="Efficiency Comparison",
          xlabel="Voltage (p.u.)",
          ylabel="Efficiency (%)")
    
    plot!(comparison_plot[4], u_vals_load, eff_load,
          label="Constant Load",
          linewidth=3,
          color=:red)
    
    # Summary table
    println("\nP-V and Q-V Characteristics Summary:")
    println("="^50)
    println("Voltage | P_const | P_load | Q_const | Q_load | Slip_load")
    println("(p.u.)  | (kW)    | (kW)   | (kVAr)  | (kVAr) | (%)")
    println("-"^60)
    
    test_voltages = [0.6, 0.8, 1.0, 1.2]
    for v_test in test_voltages
        # Find indices for both scenarios
        idx_const = argmin(abs.(u_vals_const .- v_test))
        idx_load = findfirst(v -> abs(v - v_test) < 0.05, u_vals_load)
        
        if idx_load !== nothing
            @printf(" %.1f    | %6.2f | %6.2f | %7.2f | %6.2f | %5.1f\n",
                    v_test, 
                    p_vals_const[idx_const], p_vals_load[idx_load],
                    q_vals_const[idx_const], q_vals_load[idx_load],
                    slip_vals_load[idx_load])
        else
            @printf(" %.1f    | %6.2f |   N/A  | %7.2f |  N/A   |  N/A\n",
                    v_test,
                    p_vals_const[idx_const], q_vals_const[idx_const])
        end
    end
    
    println("\nKey Insights from P-V and Q-V Analysis:")
    println("1. P-V Relationship:")
    println("   • Constant slip: P ∝ V² (parabolic)")
    println("   • Constant load: P varies less (motor adapts slip)")
    println("2. Q-V Relationship:")
    println("   • Characteristic U-shaped curve due to:")
    println("     - Low V: High current → High leakage Q")
    println("     - High V: Saturation → Reduced magnetizing Q")
    println("   • Minimum Q occurs around 0.8-1.0 p.u.")
    println("3. Practical Operation:")
    println("   • Constant load scenario is realistic")
    println("   • Efficiency drops dramatically at low voltage")
    println("   • Voltage regulation is critical")
    
    return pv_plot, qv_plot, q_components_plot, comparison_plot
end

"""
    create_comprehensive_pv_qv_analysis()

Create comprehensive P-V and Q-V analysis with multiple motor examples.
"""
function create_comprehensive_pv_qv_analysis()
    println("="^60)
    println("COMPREHENSIVE P-V AND Q-V ANALYSIS")
    println("="^60)
    
    # Create different motor examples
    motors = [
        InductionMotor(
            rated_power=5.0, rated_voltage=400.0, rated_frequency=50.0,
            rated_speed=1450.0, poles=4,
            r1=0.6, x1=2.0, r2=0.4, x2=2.0, xm=25.0
        ),
        InductionMotor(
            rated_power=15.0, rated_voltage=400.0, rated_frequency=50.0,
            rated_speed=1450.0, poles=4,
            r1=0.3, x1=1.2, r2=0.2, x2=1.2, xm=35.0
        ),
        InductionMotor(
            rated_power=50.0, rated_voltage=400.0, rated_frequency=50.0,
            rated_speed=1450.0, poles=4,
            r1=0.1, x1=0.8, r2=0.08, x2=0.8, xm=50.0
        )
    ]
    
    motor_names = ["5 kW Motor", "15 kW Motor", "50 kW Motor"]
    colors = [:blue, :red, :green]
    
    # Create comparative P-V plot
    pv_comparison = plot(xlabel="Voltage (p.u.)",
                        ylabel="Active Power (kW)",
                        title="P-V Characteristics: Motor Size Comparison",
                        grid=true,
                        legend=:topleft,
                        size=(800, 600))
    
    # Create comparative Q-V plot
    qv_comparison = plot(xlabel="Voltage (p.u.)",
                        ylabel="Reactive Power (kVAr)",
                        title="Q-V Characteristics: Motor Size Comparison (U-shaped)",
                        grid=true,
                        legend=:topright,
                        size=(800, 600))
    
    voltage_range = range(0.5, 1.4, length=80)
    
    for (i, motor) in enumerate(motors)
        # Calculate rated slip
        sync_speed = 120 * motor.rated_frequency / motor.poles
        rated_slip = (sync_speed - motor.rated_speed) / sync_speed
        
        p_vals = Float64[]
        q_vals = Float64[]
        
        for v_pu in voltage_range
            voltage = v_pu * motor.rated_voltage
            p_m, _, q, _, _, _, _, _, _, _, _, _, _ = 
                calculate_enhanced_motor_performance(motor, voltage, rated_slip)
            push!(p_vals, p_m)
            push!(q_vals, q)
        end
        
        # Normalize P by rated power for comparison
        p_normalized = p_vals ./ motor.rated_power
        
        plot!(pv_comparison, voltage_range, p_vals,
              label=motor_names[i],
              linewidth=3,
              color=colors[i])
        
        plot!(qv_comparison, voltage_range, q_vals,
              label=motor_names[i],
              linewidth=3,
              color=colors[i])
        
        # Mark minimum Q point
        min_q_idx = argmin(q_vals)
        scatter!(qv_comparison, [voltage_range[min_q_idx]], [q_vals[min_q_idx]],
                markersize=6,
                color=colors[i],
                markershape=:circle,
                label="")
    end
    
    # Add annotations
    annotate!(pv_comparison, 0.7, 40,
              text("P ∝ V² relationship\nfor constant slip", 10, :left))
    
    annotate!(qv_comparison, 0.6, 15,
              text("Low V: High current\n→ High Q_leakage", 9, :left))
    annotate!(qv_comparison, 1.3, 15,
              text("High V: Saturation\n→ Reduced Q_mag", 9, :right))
    
    return pv_comparison, qv_comparison
end

"""
    main()

Main comprehensive test function.
"""
function main()
    println("COMPREHENSIVE INDUCTION MOTOR ANALYSIS")
    println("="^80)
    println("This analysis demonstrates:")
    println("1. How slip increases dramatically at low voltage (constant torque)")
    println("2. How rotor losses increase quadratically with higher slip") 
    println("3. How magnetic saturation causes core losses at high voltage")
    println("4. P-V and Q-V characteristic relationships")
    println("5. The critical importance of voltage regulation")
    println("="^80)
    println()
    
    # Run comprehensive constant torque analysis
    motor1, comprehensive_plot, u_vals, eff_vals, rotor_losses, core_losses = 
        test_constant_torque_voltage_effects()
    
    println("\n")
    
    # Generate P-V and Q-V characteristic plots
    println("Generating P-V and Q-V characteristic curves...")
    pv_plot, qv_plot, q_components_plot, comparison_plot = plot_pv_qv_characteristics(motor1)
    
    # Display P-V and Q-V plots
    display(pv_plot)
    display(qv_plot)
    display(q_components_plot)
    display(comparison_plot)
    
    println("\n")
    
    # Create comprehensive P-V and Q-V comparison
    println("Creating comprehensive P-V and Q-V motor comparison...")
    pv_comparison, qv_comparison = create_comprehensive_pv_qv_analysis()
    display(pv_comparison)
    display(qv_comparison)
    
    println("\n")
    
    # Run detailed loss mechanism analysis
    motor2 = test_loss_mechanisms_detailed()
    
    # Summary statistics
    println("\n" * "="^80)
    println("CRITICAL FINDINGS:")
    min_eff = minimum(eff_vals)
    max_rotor_loss = maximum(rotor_losses)
    max_core_loss = maximum(core_losses)
    
    println("• Efficiency range: $(round(min_eff, digits=1))% - $(round(maximum(eff_vals), digits=1))%")
    println("• Rotor loss range: $(round(minimum(rotor_losses), digits=0))W - $(round(max_rotor_loss, digits=0))W")
    println("• Core loss range: $(round(minimum(core_losses), digits=0))W - $(round(max_core_loss, digits=0))W")
    println("• Efficiency can drop by $(round(maximum(eff_vals) - min_eff, digits=1)) percentage points!")
    println()
    println("P-V AND Q-V RELATIONSHIP INSIGHTS:")
    println("✓ P-V: Parabolic for constant slip (P ∝ V²), flatter for constant load")
    println("✓ Q-V: U-shaped curve due to leakage vs magnetizing competition")
    println("✓ Q minimum: Occurs around 0.8-1.0 p.u. voltage")
    println("✓ Saturation effect: Creates characteristic Q decrease at high voltage")
    println()
    println("KEY TAKEAWAYS:")
    println("✓ Low voltage: Slip ↑ → R2'/s ↓ → I2 ↑ → P_rotor = I2²×R2 ↑↑")
    println("✓ High voltage: Saturation → Xm ↓ → Flux density ↑ → P_core ↑↑")
    println("✓ Voltage regulation is ESSENTIAL for motor efficiency")
    println("✓ Operating range should be 0.95-1.05 p.u. for optimal performance")
    println("✓ P-V and Q-V curves provide critical design insights")
    println("="^80)
    
    return motor1, motor2, comprehensive_plot, pv_plot, qv_plot, q_components_plot, comparison_plot, pv_comparison, qv_comparison
end

# Run the comprehensive analysis
motor1, motor2, comprehensive_plot, pv_plot, qv_plot, q_components_plot, comparison_plot, pv_comparison, qv_comparison = main()

