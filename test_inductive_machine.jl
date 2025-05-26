include("visualization/induction_motor_analysis.jl")

"""
    test_upq_relationship()

Test and visualize the relationship between voltage (U), active power (P), 
and reactive power (Q) for induction motors under different operating conditions.
"""
function test_upq_relationship()
    println("="^60)
    println("TESTING U-P-Q RELATIONSHIP FOR INDUCTION MOTOR")
    println("="^60)
    
    # Create a test motor
    motor = InductionMotor(
        rated_power=10.0,      # 10 kW
        rated_voltage=400.0,   # 400 V
        rated_frequency=50.0,  # 50 Hz
        rated_speed=1450.0,    # 1450 rpm
        poles=4,
        r1=0.5,               # Stator resistance
        x1=1.0,               # Stator reactance
        r2=0.3,               # Rotor resistance
        x2=1.0,               # Rotor reactance
        xm=30.0               # Magnetizing reactance
    )
    
    # Calculate synchronous speed and rated slip
    sync_speed = 120 * motor.rated_frequency / motor.poles
    rated_slip = (sync_speed - motor.rated_speed) / sync_speed
    
    println("Motor Parameters:")
    println("  Rated Power: $(motor.rated_power) kW")
    println("  Rated Voltage: $(motor.rated_voltage) V")
    println("  Rated Slip: $(round(rated_slip*100, digits=2))%")
    println()
    
    # Test 1: U-P-Q relationship at constant slip (rated slip)
    println("Test 1: U-P-Q at Constant Slip ($(round(rated_slip*100, digits=1))%)")
    println("-"^50)
    
    voltage_range = range(0.7, 1.3, length=50)
    u_vals = Float64[]
    p_vals = Float64[]
    q_vals = Float64[]
    
    for v_pu in voltage_range
        voltage = v_pu * motor.rated_voltage
        p_m, _, q, _, _, _ = calculate_motor_performance(motor, voltage, rated_slip)
        push!(u_vals, v_pu)
        push!(p_vals, p_m)
        push!(q_vals, q)
    end
    
    # Plot U-P and U-Q relationships
    p1 = plot(u_vals, p_vals, 
              label="P(U)", 
              xlabel="Voltage (p.u.)", 
              ylabel="Power (kW)",
              title="Active Power vs Voltage",
              linewidth=2,
              color=:blue,
              grid=true)
    
    p2 = plot(u_vals, q_vals,
              label="Q(U)",
              xlabel="Voltage (p.u.)",
              ylabel="Reactive Power (kVAr)",
              title="Reactive Power vs Voltage", 
              linewidth=2,
              color=:red,
              grid=true)
    
    # P-Q curve at constant slip
    p3 = plot(p_vals, q_vals,
              label="Constant Slip Curve",
              xlabel="Active Power (kW)",
              ylabel="Reactive Power (kVAr)",
              title="P-Q Relationship at Constant Slip",
              linewidth=2,
              color=:green,
              grid=true)
    
    # Mark rated point
    p_rated, _, q_rated, _, _, _ = calculate_motor_performance(motor, motor.rated_voltage, rated_slip)
    scatter!(p3, [p_rated], [q_rated], 
            markersize=8, 
            color=:red, 
            label="Rated Point",
            markershape=:star)
    
    # Test 2: U-P-Q relationship at constant power
    println("Test 2: U-Q at Constant Power")
    println("-"^50)
    
    target_powers = [0.3, 0.5, 0.7, 1.0] .* motor.rated_power
    p4 = plot(xlabel="Voltage (p.u.)",
              ylabel="Reactive Power (kVAr)",
              title="Q vs U at Different Constant Powers",
              legend=:topright,
              grid=true)
    
    colors = [:blue, :red, :green, :purple]
    
    for (i, target_p) in enumerate(target_powers)
        u_const_p = Float64[]
        q_const_p = Float64[]
        
        for v_pu in voltage_range
            voltage = v_pu * motor.rated_voltage
            
            # Find slip that gives target power at this voltage
            slip_found = false
            for s in range(0.001, 0.5, length=100)
                p_m, _, q, _, _, _ = calculate_motor_performance(motor, voltage, s)
                if abs(p_m - target_p) < 0.1  # Within 0.1 kW tolerance
                    push!(u_const_p, v_pu)
                    push!(q_const_p, q)
                    slip_found = true
                    break
                end
            end
        end
        
        if !isempty(u_const_p)
            plot!(p4, u_const_p, q_const_p,
                  label="P = $(round(target_p, digits=1)) kW",
                  linewidth=2,
                  color=colors[i])
        end
    end
    
    # Test 3: 3D visualization of U-P-Q relationship
    println("Test 3: 3D U-P-Q Surface")
    println("-"^50)
    
    # Create 3D surface data
    voltage_3d = range(0.7, 1.3, length=20)
    slip_3d = range(0.01, 0.4, length=20)
    
    U_surface = Float64[]
    P_surface = Float64[]
    Q_surface = Float64[]
    
    for v_pu in voltage_3d
        for s in slip_3d
            voltage = v_pu * motor.rated_voltage
            p_m, _, q, _, _, _ = calculate_motor_performance(motor, voltage, s)
            push!(U_surface, v_pu)
            push!(P_surface, p_m)
            push!(Q_surface, q)
        end
    end
    
    p5 = scatter(U_surface, P_surface, Q_surface,
                 xlabel="Voltage (p.u.)",
                 ylabel="Active Power (kW)",
                 zlabel="Reactive Power (kVAr)",
                 title="3D U-P-Q Relationship",
                 markersize=2,
                 markercolor=:viridis,
                 camera=(45, 30),
                 size=(700, 600))
    
    # Mark rated operating point in 3D
    scatter!(p5, [1.0], [p_rated], [q_rated],
            markersize=8,
            color=:red,
            markershape=:star,
            label="Rated Point")
    
    # Test 4: Contour plot of U-P-Q relationship
    println("Test 4: U-P Contour Plot with Q Lines")
    println("-"^50)
    
    # Create grid for contour plot
    u_grid = range(0.7, 1.3, length=30)
    s_grid = range(0.01, 0.4, length=30)
    
    P_grid = zeros(length(s_grid), length(u_grid))
    Q_grid = zeros(length(s_grid), length(u_grid))
    
    for (i, s) in enumerate(s_grid)
        for (j, v_pu) in enumerate(u_grid)
            voltage = v_pu * motor.rated_voltage
            p_m, _, q, _, _, _ = calculate_motor_performance(motor, voltage, s)
            P_grid[i, j] = p_m
            Q_grid[i, j] = q
        end
    end
    
    p6 = contour(u_grid, s_grid, P_grid,
                xlabel="Voltage (p.u.)",
                ylabel="Slip",
                title="Power Contours with Q Lines",
                fill=true,
                color=:viridis,
                size=(600, 500))
    
    # Add Q contour lines
    contour!(p6, u_grid, s_grid, Q_grid,
            levels=10,
            color=:black,
            linewidth=1,
            linestyle=:dash,
            label="Q contours")
    
    # Display all plots
    combined_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
    
    display(combined_plot)
    display(p5)
    display(p6)
    
    # Print some numerical results
    println("\nNumerical Results at Different Voltages (Rated Slip):")
    println("Voltage(p.u.) | Power(kW) | Q(kVAr) | P/Q Ratio")
    println("-"^50)
    for v_pu in [0.8, 0.9, 1.0, 1.1, 1.2]
        voltage = v_pu * motor.rated_voltage
        p_m, _, q, _, _, _ = calculate_motor_performance(motor, voltage, rated_slip)
        pq_ratio = q != 0 ? p_m/q : Inf
        @printf("    %.1f     |   %.2f   |  %.2f  |   %.2f\n", v_pu, p_m, q, pq_ratio)
    end
    
    println("\n" * "="^60)
    println("U-P-Q RELATIONSHIP TESTING COMPLETED")
    println("="^60)
    
    return motor, combined_plot, p5, p6
end

# Run the test directly
motor, combined_plot, surface_plot, contour_plot = test_upq_relationship()
