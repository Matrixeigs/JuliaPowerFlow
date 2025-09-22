using DifferentialEquations
using Plots

# Generator parameters (with some corrections)
struct GeneratorParams
    H::Float64      # Inertia constant (s)
    D::Float64      # Damping coefficient (pu)
    ωs::Float64     # Synchronous speed (rad/s)
    
    # Reactances (pu) - slightly adjusted to avoid numerical issues
    Xd::Float64     # d-axis synchronous reactance
    Xq::Float64     # q-axis synchronous reactance
    X_d_prime::Float64    # d-axis transient reactance
    X_q_prime::Float64    # q-axis transient reactance
    X_d_double::Float64   # d-axis subtransient reactance
    X_q_double::Float64   # q-axis subtransient reactance (made slightly different)
    Xl::Float64     # Leakage reactance
    
    # Time constants (s)
    T_d0_prime::Float64   # d-axis transient time constant
    T_q0_prime::Float64   # q-axis transient time constant
    T_d0_double::Float64  # d-axis subtransient time constant
    T_q0_double::Float64  # q-axis subtransient time constant
    T_d_double::Float64   # d-axis subtransient time constant (loaded)
    T_q_double::Float64   # q-axis subtransient time constant (loaded)
    
    # Network parameters
    Vt::Float64     # Terminal voltage magnitude (pu)
    Xe::Float64     # External reactance (pu)
    Pe::Float64     # Electrical power output (pu)
end

# Initialize parameters with fixes
function get_params()
    return GeneratorParams(
        3.5,    # H
        2.0,    # D
        377.0,  # ωs
        1.8,    # Xd
        1.7,    # Xq
        0.3,    # X_d_prime
        0.55,   # X_q_prime
        0.25,   # X_d_double
        0.24,   # X_q_double (made slightly different to avoid algebraic issues)
        0.2,    # Xl
        8.0,    # T_d0_prime
        0.4,    # T_q0_prime
        0.03,   # T_d0_double
        0.05,   # T_q0_double
        0.035,  # T_d_double
        0.06,   # T_q_double
        1.0,    # Vt
        0.4,    # Xe
        0.9     # Pe
    )
end

# Calculate initial conditions with error checking
function calculate_initial_conditions(p::GeneratorParams)
    # Check if power equation is solvable
    arcsin_arg = p.Pe * p.Xe / p.Vt
    if abs(arcsin_arg) > 1.0
        error("Power equation not solvable: Pe*Xe/Vt = $arcsin_arg > 1.0")
    end
    
    δ0 = asin(arcsin_arg)
    ω0 = 0.0  # Initial frequency deviation
    
    # Calculate initial flux linkages with numerical conditioning
    cos_δ = cos(δ0)
    sin_δ = sin(δ0)
    
    # Add small epsilon to prevent division by zero
    ε = 1e-12
    
    # Initial conditions for flux linkages (simplified steady-state)
    ψ_d_prime0 = p.Vt * cos_δ
    ψ_q_prime0 = p.Vt * sin_δ
    ψ_d_double0 = p.Vt * cos_δ
    ψ_q_double0 = p.Vt * sin_δ
    
    return [δ0, ω0, ψ_d_prime0, ψ_q_prime0, ψ_d_double0, ψ_q_double0]
end

# Differential equation system with numerical conditioning
function generator_dynamics!(du, u, p, t)
    δ, ω, ψ_d_prime, ψ_q_prime, ψ_d_double, ψ_q_double = u
    
    # Add small epsilon for numerical stability
    ε = 1e-12
    
    # Calculate currents with conditioning
    cos_δ = cos(δ)
    sin_δ = sin(δ)
    
    # Current calculations with numerical conditioning
    denom_d = (p.X_d_double + p.Xe) + ε
    denom_q = (p.X_q_double + p.Xe) + ε
    
    i_d = (p.Vt * cos_δ - ψ_d_double) / denom_d
    i_q = (p.Vt * sin_δ - ψ_q_double) / denom_q
    
    # Electrical power
    Pe = ψ_d_double * i_q - ψ_q_double * i_d
    
    # Mechanical equation
    du[1] = p.ωs * ω  # dδ/dt
    du[2] = (1/(2*p.H)) * (p.Pe - Pe - p.D * ω)  # dω/dt
    
    # Flux linkage equations with conditioning
    du[3] = (1/p.T_d0_prime) * (-ψ_d_prime + (p.Xd - p.X_d_prime) * i_d)
    du[4] = (1/p.T_q0_prime) * (-ψ_q_prime + (p.Xq - p.X_q_prime) * i_q)
    
    du[5] = (1/p.T_d_double) * (-ψ_d_double + ψ_d_prime - (p.X_d_prime - p.X_d_double) * i_d)
    du[6] = (1/p.T_q_double) * (-ψ_q_double + ψ_q_prime - (p.X_q_prime - p.X_q_double) * i_q)
    
    return nothing
end

# Main simulation function
function simulate_generator()
    p = get_params()
    
    println("Generator Parameters:")
    println("Stiffness ratio: ", p.T_d0_prime / min(p.T_d0_double, p.T_q0_double))
    
    # Calculate initial conditions
    u0 = calculate_initial_conditions(p)
    println("Initial conditions: ", u0)
    
    # Time span
    tspan = (0.0, 2.0)  # Shorter simulation time initially
    
    # Create ODE problem
    prob = ODEProblem(generator_dynamics!, u0, tspan, p)
    
    # Use stiff solver with appropriate tolerances
    sol = solve(prob, Rodas4(), 
                reltol=1e-6, abstol=1e-8,
                maxiters=1e6,
                dtmax=0.001)  # Limit maximum time step
    
    return sol, p
end

# Run simulation and plot results
try
    sol, p = simulate_generator()
    
    println("Simulation completed successfully!")
    println("Solution length: ", length(sol.t))
    println("Final time: ", sol.t[end])
    
    # Plot results
    plot_δ = plot(sol.t, sol[1,:] * 180/π, 
                  title="Rotor Angle", xlabel="Time (s)", ylabel="δ (degrees)")
    
    plot_ω = plot(sol.t, sol[2,:], 
                  title="Frequency Deviation", xlabel="Time (s)", ylabel="ω (pu)")
    
    plot_flux = plot(sol.t, [sol[3,:] sol[4,:] sol[5,:] sol[6,:]], 
                     title="Flux Linkages", xlabel="Time (s)", ylabel="Flux (pu)",
                     label=["ψ'd" "ψ'q" "ψ''d" "ψ''q"])
    
    plot(plot_δ, plot_ω, plot_flux, layout=(3,1), size=(800,600))
    
catch e
    println("Error in simulation: ", e)
    println("Try the following fixes:")
    println("1. Reduce time step further")
    println("2. Check parameter values")
    println("3. Use different initial conditions")
    println("4. Try different solver (e.g., Rosenbrock23(), TRBDF2())")
end

# Alternative solver options to try if Rodas4() doesn't work:
# solve(prob, Rosenbrock23(), reltol=1e-6, abstol=1e-8, dtmax=0.001)
# solve(prob, TRBDF2(), reltol=1e-6, abstol=1e-8, dtmax=0.001)
# solve(prob, RadauIIA5(), reltol=1e-6, abstol=1e-8, dtmax=0.001)
# solve(prob, KenCarp4(), reltol=1e-6, abstol=1e-8, dtmax=0.001)