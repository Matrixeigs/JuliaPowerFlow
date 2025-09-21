using DifferentialEquations
using Plots
using LinearAlgebra

# Synchronous Generator Seventh-Order Model
# State variables: [δ, ω, E'q, E'd, E''q, E''d, ψ2q]
# where:
# δ - rotor angle (rad)
# ω - rotor speed deviation (pu)
# E'q - q-axis transient EMF (pu)
# E'd - d-axis transient EMF (pu)
# E''q - q-axis subtransient EMF (pu)
# E''d - d-axis subtransient EMF (pu)
# ψ2q - q-axis damper winding flux linkage (pu)

struct GeneratorParameters
    # Machine parameters (per unit)
    H::Float64          # Inertia constant (s)
    D::Float64          # Damping coefficient (pu)
    ωs::Float64         # Synchronous speed (rad/s)
    
    # Reactances (pu)
    Xd::Float64         # d-axis synchronous reactance
    Xq::Float64         # q-axis synchronous reactance
    X_d_prime::Float64  # d-axis transient reactance
    X_q_prime::Float64  # q-axis transient reactance
    X_d_double::Float64 # d-axis subtransient reactance
    X_q_double::Float64 # q-axis subtransient reactance
    Xl::Float64         # Leakage reactance
    
    # Time constants (s)
    T_d0_prime::Float64 # d-axis transient time constant
    T_q0_prime::Float64 # q-axis transient time constant
    T_d0_double::Float64# d-axis subtransient time constant
    T_q0_double::Float64# q-axis subtransient time constant
    T_d_double::Float64 # d-axis subtransient time constant (loaded)
    T_q_double::Float64 # q-axis subtransient time constant (loaded)
    
    # Network parameters
    Vt::Float64         # Terminal voltage magnitude (pu)
    Xe::Float64         # External reactance (pu)
    Pe::Float64         # Electrical power output (pu)
end

# Default generator parameters (typical values) - Fixed constructor
function default_generator_params()
    return GeneratorParameters(
        3.5,            # H - Inertia constant
        2.0,            # D - Damping coefficient
        377.0,          # ωs - 60 Hz synchronous speed
        
        1.8,            # Xd - d-axis synchronous reactance
        1.7,            # Xq - q-axis synchronous reactance
        0.3,            # X_d_prime - d-axis transient reactance
        0.55,           # X_q_prime - q-axis transient reactance
        0.25,           # X_d_double - d-axis subtransient reactance
        0.25,           # X_q_double - q-axis subtransient reactance
        0.2,            # Xl - Leakage reactance
        
        8.0,            # T_d0_prime - d-axis transient time constant
        0.4,            # T_q0_prime - q-axis transient time constant
        0.03,           # T_d0_double - d-axis subtransient time constant
        0.05,           # T_q0_double - q-axis subtransient time constant
        0.035,          # T_d_double - d-axis subtransient time constant (loaded)
        0.06,           # T_q_double - q-axis subtransient time constant (loaded)
        
        1.0,            # Vt - Terminal voltage
        0.4,            # Xe - External reactance
        0.9             # Pe - Electrical power output
    )
end

# Calculate steady-state initial conditions
function calculate_initial_conditions(params::GeneratorParameters)
    # Calculate initial rotor angle and currents
    δ0 = asin(params.Pe * params.Xe / params.Vt)
    
    # Terminal voltage components
    Vd = params.Vt * sin(δ0)
    Vq = params.Vt * cos(δ0)
    
    # Current components (simplified)
    Id = (params.Pe * params.Xe * cos(δ0)) / (params.Vt^2)
    Iq = (params.Pe * params.Xe * sin(δ0)) / (params.Vt^2)
    
    # Initial EMF values
    Eq_prime0 = Vq + params.X_d_prime * Id
    Ed_prime0 = Vd - params.X_q_prime * Iq
    Eq_double0 = Vq + params.X_d_double * Id
    Ed_double0 = Vd - params.X_q_double * Iq
    
    # Field flux linkage
    ψ2q0 = Ed_prime0
    
    return [δ0, 0.0, Eq_prime0, Ed_prime0, Eq_double0, Ed_double0, ψ2q0]
end

# Seventh-order differential equation system
function generator_dynamics!(dx, x, p, t)
    params, Pm, Efd = p
    δ, ω, Eq_prime, Ed_prime, Eq_double, Ed_double, ψ2q = x
    
    # Terminal voltage components in dq frame
    Vd = params.Vt * sin(δ)
    Vq = params.Vt * cos(δ)
    
    # Calculate currents using subtransient model
    # Network equations with external impedance
    Xt = params.Xe + params.X_d_double
    Yt = params.Xe + params.X_q_double
    
    denom = Xt^2 + Yt^2
    Id = (Xt * (Eq_double - Vq) + Yt * (Ed_double - Vd)) / denom
    Iq = (Yt * (Eq_double - Vq) - Xt * (Ed_double - Vd)) / denom
    
    # Electrical power
    Pe = Vd * Id + Vq * Iq
    
    # Mechanical equation (swing equation)
    dx[1] = params.ωs * ω  # dδ/dt
    dx[2] = (Pm - Pe - params.D * ω) / (2 * params.H)  # dω/dt
    
    # Flux linkage dynamics
    # d-axis transient EMF
    dx[3] = (Efd - Eq_prime + (params.Xd - params.X_d_prime) * Id) / params.T_d0_prime
    
    # q-axis transient EMF  
    dx[4] = (-Ed_prime - (params.Xq - params.X_q_prime) * Iq) / params.T_q0_prime
    
    # d-axis subtransient EMF
    dx[5] = (Eq_prime - Eq_double + (params.X_d_prime - params.X_d_double) * Id) / params.T_d_double
    
    # q-axis subtransient EMF
    dx[6] = (Ed_prime - Ed_double - (params.X_q_prime - params.X_q_double) * Iq) / params.T_q_double
    
    # q-axis damper winding flux linkage
    dx[7] = (-ψ2q + Ed_prime) / params.T_q0_double
    
    return nothing
end

# Simulation function
function simulate_generator(; 
    params = default_generator_params(),
    Pm = 0.9,           # Mechanical power input (pu)
    Efd = 1.2,          # Field voltage (pu)
    disturbance_time = 1.0,  # Time of disturbance (s)
    disturbance_magnitude = 0.1,  # Magnitude of power disturbance
    tspan = (0.0, 5.0)  # Simulation time span
    )
    
    # Calculate initial conditions
    x0 = calculate_initial_conditions(params)
    
    # Define time-varying inputs (mechanical power disturbance)
    function inputs(t)
        if t >= disturbance_time && t <= disturbance_time + 0.1
            return Pm + disturbance_magnitude, Efd
        else
            return Pm, Efd
        end
    end
    
    # Create parameter tuple
    p = (params, Pm, Efd)
    
    # Solve ODE with time-varying parameters
    function generator_with_disturbance!(dx, x, p, t)
        params, _, _ = p
        Pm_t, Efd_t = inputs(t)
        p_new = (params, Pm_t, Efd_t)
        generator_dynamics!(dx, x, p_new, t)
    end
    
    # Define and solve the problem
    prob = ODEProblem(generator_with_disturbance!, x0, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-10)
    
    return sol, params
end

# Post-processing function to calculate additional quantities
function calculate_outputs(sol, params)
    t = sol.t
    n = length(t)
    
    # Initialize output arrays
    Pe = zeros(n)
    Id = zeros(n)
    Iq = zeros(n)
    Vd = zeros(n)
    Vq = zeros(n)
    Te = zeros(n)
    
    for i in 1:n
        δ, ω, Eq_prime, Ed_prime, Eq_double, Ed_double, ψ2q = sol.u[i]
        
        # Terminal voltage components
        Vd[i] = params.Vt * sin(δ)
        Vq[i] = params.Vt * cos(δ)
        
        # Calculate currents
        Xt = params.Xe + params.X_d_double
        Yt = params.Xe + params.X_q_double
        denom = Xt^2 + Yt^2
        
        Id[i] = (Xt * (Eq_double - Vq[i]) + Yt * (Ed_double - Vd[i])) / denom
        Iq[i] = (Yt * (Eq_double - Vq[i]) - Xt * (Ed_double - Vd[i])) / denom
        
        # Electrical power and torque
        Pe[i] = Vd[i] * Id[i] + Vq[i] * Iq[i]
        Te[i] = Eq_double * Iq[i] + Ed_double * Id[i]
    end
    
    return (t=t, Pe=Pe, Id=Id, Iq=Iq, Vd=Vd, Vq=Vq, Te=Te)
end

# Main simulation and plotting
function main_simulation()
    println("Starting Synchronous Generator Simulation...")
    
    # Run simulation
    sol, params = simulate_generator()
    outputs = calculate_outputs(sol, params)
    
    # Extract state variables
    δ = [u[1] for u in sol.u] * 180/π  # Convert to degrees
    ω = [u[2] for u in sol.u]
    Eq_prime = [u[3] for u in sol.u]
    Ed_prime = [u[4] for u in sol.u]
    Eq_double = [u[5] for u in sol.u]
    Ed_double = [u[6] for u in sol.u]
    
    # Create comprehensive plots
    p1 = plot(sol.t, δ, title="Rotor Angle", xlabel="Time (s)", ylabel="δ (degrees)", 
              linewidth=2, color=:blue)
    
    p2 = plot(sol.t, ω, title="Speed Deviation", xlabel="Time (s)", ylabel="Δω (pu)", 
              linewidth=2, color=:red)
    
    p3 = plot(sol.t, [Eq_prime, Ed_prime], title="Transient EMFs", 
              xlabel="Time (s)", ylabel="EMF (pu)", 
              label=["E'q" "E'd"], linewidth=2)
    
    p4 = plot(sol.t, [Eq_double, Ed_double], title="Subtransient EMFs", 
              xlabel="Time (s)", ylabel="EMF (pu)", 
              label=["E''q" "E''d"], linewidth=2)
    
    p5 = plot(outputs.t, outputs.Pe, title="Electrical Power", 
              xlabel="Time (s)", ylabel="Pe (pu)", linewidth=2, color=:green)
    
    p6 = plot(outputs.t, [outputs.Id, outputs.Iq], title="Stator Currents", 
              xlabel="Time (s)", ylabel="Current (pu)", 
              label=["Id" "Iq"], linewidth=2)
    
    # Combine all plots
    combined_plot = plot(p1, p2, p3, p4, p5, p6, layout=(3,2), size=(1000, 800))
    display(combined_plot)
    
    println("Simulation completed successfully!")
    println("Initial rotor angle: $(round(δ[1], digits=2)) degrees")
    println("Steady-state power: $(round(outputs.Pe[end], digits=3)) pu")
    
    return sol, outputs, params
end

# Additional analysis functions
function frequency_analysis(sol)
    ω = [u[2] for u in sol.u]
    freq = 60 .+ ω .* 60  # Convert to Hz
    
    freq_plot = plot(sol.t, freq, title="System Frequency", 
         xlabel="Time (s)", ylabel="Frequency (Hz)", 
         linewidth=2, color=:purple)
    display(freq_plot)
    return freq_plot
end

function electromagnetic_analysis(sol, params)
    # Calculate electromagnetic quantities
    n = length(sol.t)
    flux_d = zeros(n)
    flux_q = zeros(n)
    field_current = zeros(n)
    
    for i in 1:n
        δ, ω, Eq_prime, Ed_prime, Eq_double, Ed_double, ψ2q = sol.u[i]
        
        # Flux linkages
        flux_d[i] = Ed_double
        flux_q[i] = Eq_double
        
        # Approximate field current (simplified)
        field_current[i] = Eq_prime / params.Xd
    end
    
    p1 = plot(sol.t, [flux_d, flux_q], title="Flux Linkages", 
              xlabel="Time (s)", ylabel="Flux (pu)", 
              label=["ψd" "ψq"], linewidth=2)
    
    p2 = plot(sol.t, field_current, title="Field Current", 
              xlabel="Time (s)", ylabel="If (pu)", 
              linewidth=2, color=:orange)
    
    combined = plot(p1, p2, layout=(2,1), size=(800, 600))
    display(combined)
    return combined
end

# Run the main simulation
println("=== Synchronous Generator Seventh-Order Model Simulation ===")
sol, outputs, params = main_simulation()

# Additional analyses
println("\n=== Frequency Analysis ===")
freq_plot = frequency_analysis(sol)

println("\n=== Electromagnetic Analysis ===")
em_plot = electromagnetic_analysis(sol, params)

println("\n=== Generator Parameters Summary ===")
println("Inertia constant H: $(params.H) s")
println("Damping coefficient D: $(params.D) pu")
println("d-axis reactances: Xd=$(params.Xd), X'd=$(params.X_d_prime), X''d=$(params.X_d_double)")
println("q-axis reactances: Xq=$(params.Xq), X'q=$(params.X_q_prime), X''q=$(params.X_q_double)")
println("Time constants: T'd0=$(params.T_d0_prime)s, T'q0=$(params.T_q0_prime)s")
println("Subtransient time constants: T''d0=$(params.T_d0_double)s, T''q0=$(params.T_q0_double)s")

println("\n=== Simulation Results Summary ===")
println("Maximum rotor angle swing: $(round(maximum(abs.([u[1] for u in sol.u])) * 180/π, digits=2)) degrees")
println("Maximum speed deviation: $(round(maximum(abs.([u[2] for u in sol.u])), digits=4)) pu")
println("Settling time (approx): $(round(sol.t[end], digits=1)) seconds")