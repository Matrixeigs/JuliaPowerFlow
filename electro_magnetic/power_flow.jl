# power_flow.jl - Power Flow Solver for System Initialization

using LinearAlgebra, SparseArrays

"""
    run_power_flow(components)

Runs a Newton-Raphson power flow to find the steady-state solution.
"""
function run_power_flow(components)
    println("  - Running Newton-Raphson power flow...")
    
    max_iter = 20
    tolerance = 1e-6
    
    # Extract parameters
    gen = components.generator
    pv = components.pv_system
    bess = components.bess
    motor = components.motor
    ev_charger = components.ev_charger
    
    # Initial guess for voltages (flat start)
    V = ones(ComplexF64, NUM_BUSES)
    
    for iter in 1:max_iter
        # Calculate power injections
        S_gen = gen.P_out + 1im * gen.Q_out
        S_pv = pv.P_out + 1im * pv.Q_out
        S_bess = bess.P_out + 1im * bess.Q_out
        
        S_load_motor = motor.P_in + 1im * motor.Q_in
        S_load_ev = ev_charger.P_in + 1im * ev_charger.Q_in
        
        S_injected = zeros(ComplexF64, NUM_BUSES)
        S_injected[1] = S_gen  # Slack bus
        S_injected[2] = S_pv + S_bess
        S_injected[3] = -S_load_motor - S_load_ev
        
        # Build Admittance Matrix (Ybus)
        Ybus = zeros(ComplexF64, NUM_BUSES, NUM_BUSES)
        Ybus[1,1] = 1.0/components.line.Z; Ybus[1,2] = -1.0/components.line.Z;
        Ybus[2,1] = -1.0/components.line.Z; Ybus[2,2] = 1.0/components.line.Z + 1.0/components.transformer.Z;
        Ybus[2,3] = -1.0/components.transformer.Z; Ybus[3,2] = -1.0/components.transformer.Z;
        Ybus[3,3] = 1.0/components.transformer.Z;

        # Calculate power mismatch
        I = Ybus * V
        S_calc = V .* conj.(I)
        
        # Mismatch for PQ buses (bus 2 and 3)
        dS = S_injected[2:end] - S_calc[2:end]
        
        # Check for convergence
        mismatch = maximum(abs.(dS))
        if mismatch < tolerance
            println("  - Power flow converged in $iter iterations.")
            # Update component voltages based on power flow solution
            components.generator.V_t = V[1]
            components.pv_system.V_t = V[2]
            components.bess.V_t = V[2]
            components.motor.V_t = V[3]
            return V
        end
        
        # Build Jacobian (simplified for illustration)
        J = real.([Ybus[2:3, 2:3] zeros(2,2); zeros(2,2) Ybus[2:3, 2:3]])
        
        # Solve for voltage correction
        mismatch_vec = [real(dS); imag(dS)]
        dV_vec = J \ mismatch_vec
        
        # Update voltages (for PQ buses)
        V[2:end] .+= dV_vec[1:2] + 1im * dV_vec[3:4]
    end
    
    @warn "Power flow did not converge after $max_iter iterations."
    return V
end
