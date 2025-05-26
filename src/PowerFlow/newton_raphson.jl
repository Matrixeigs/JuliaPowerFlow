"""
    newton_raphson_power_flow(sys::PowerSystem; max_iter=30, tol=1e-6, verbose=true)

Solve power flow using Newton-Raphson method.
"""
function newton_raphson_power_flow(sys::PowerSystem; max_iter=30, tol=1e-6, verbose=true)
    # Build Y-bus
    Ybus, bus_to_idx = build_ybus(sys)
    
    # Initialize voltage vector
    V0 = initialize_voltage(sys, bus_to_idx)
    
    # Identify bus types
    ref, pv, pq = identify_bus_types(sys)
    ref_idx = [bus_to_idx[id] for id in ref]
    pv_idx = [bus_to_idx[id] for id in pv]
    pq_idx = [bus_to_idx[id] for id in pq]
    
    # Calculate power injection
    Sbus = calculate_power_injection(sys, bus_to_idx)
    
    # Solve using Newton-Raphson
    return newton_raphson_solve(Ybus, Sbus, V0, pv_idx, pq_idx, ref_idx, 
                               max_iter=max_iter, tol=tol, verbose=verbose)
end

"""
    newton_raphson_solve(Ybus, Sbus, V0, pv, pq, ref; max_iter=30, tol=1e-6, verbose=true)

Core Newton-Raphson solver implementation.
"""
function newton_raphson_solve(Ybus, Sbus, V0, pv, pq, ref; max_iter=30, tol=1e-6, verbose=true)
    n = length(V0)
    V = copy(V0)
    Va = angle.(V)
    Vm = abs.(V)
    iteration = 0
    converged = false
    
    # Set reference bus voltage and angle
    Va[ref] .= angle.(V0[ref])
    Vm[ref] .= abs.(V0[ref])
    
    # Set PV bus voltage magnitude
    Vm[pv] .= abs.(V0[pv])
    
    pvpq = [pv; pq]
    
    if verbose
        println("Newton-Raphson power flow started")
        println("Max iterations: $max_iter, Tolerance: $tol")
        println("Buses: $(n), PV: $(length(pv)), PQ: $(length(pq)), Ref: $(length(ref))")
    end
    
    while !converged && iteration < max_iter
        iteration += 1
        
        # Update voltage vector
        V = Vm .* exp.(im .* Va)
        
        # Calculate power injections
        Ibus = Ybus * V
        S = V .* conj.(Ibus)
        P = real.(S)
        Q = imag.(S)
        
        # Calculate mismatches
        dP = real.(Sbus) - P
        dQ = imag.(Sbus) - Q
        
        # Form mismatch vector
        F = [dP[pvpq]; dQ[pq]]
        norm_F = norm(F, Inf)
        
        if verbose
            @printf("Iteration %3d: Max mismatch = %.10f\n", iteration, norm_F)
        end
        
        if norm_F < tol
            converged = true
            if verbose
                println("Newton-Raphson converged!")
            end
            break
        end
        
        # Build Jacobian
        J = create_jacobian(Ybus, V, pv, pq)
        
        # Solve linear system
        dx = J \ F
        
        # Update state variables
        Va[pvpq] += dx[1:length(pvpq)]
        if length(pq) > 0
            Vm[pq] .*= (1.0 .+ dx[length(pvpq)+1:end])
        end
    end
    
    if !converged
        @warn "Newton-Raphson did not converge within maximum iterations!"
    end
    
    # Final voltage update
    V = Vm .* exp.(im .* Va)
    
    # Calculate final power injections
    Ibus = Ybus * V
    S = V .* conj.(Ibus)
    
    return V, converged, iteration, S
end

"""
    initialize_voltage(sys::PowerSystem, bus_to_idx::Dict)

Initialize voltage vector for the power system.
"""
function initialize_voltage(sys::PowerSystem, bus_to_idx::Dict)
    n = length(bus_to_idx)
    V = ones(Complex{Float64}, n)
    
    for (id, bus) in sys.buses
        idx = bus_to_idx[id]
        V[idx] = bus.voltage_magnitude * exp(im * bus.voltage_angle)
    end
    
    # Set generator voltages
    for (_, gen) in sys.generators
        if gen.is_online
            bus_idx = bus_to_idx[gen.bus_id]
            V[bus_idx] = gen.voltage_setpoint * exp(im * angle(V[bus_idx]))
        end
    end
    
    return V
end
