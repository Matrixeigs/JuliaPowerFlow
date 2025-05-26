"""
    build_ybus(sys::PowerSystem)

Build the bus admittance matrix for the power system.
"""
function build_ybus(sys::PowerSystem)
    nb = get_bus_count(sys)
    bus_ids = sort(collect(keys(sys.buses)))
    
    # Create bus ID to index mapping
    bus_to_idx = Dict(id => i for (i, id) in enumerate(bus_ids))
    
    # Initialize Y-bus
    Ybus = spzeros(Complex{Float64}, nb, nb)
    
    # Add branch admittances
    for (_, branch) in sys.branches
        if !branch.is_online
            continue
        end
        
        from_idx = bus_to_idx[branch.from_bus]
        to_idx = bus_to_idx[branch.to_bus]
        
        # Get branch parameters
        y = get_admittance(branch)
        b_sh = get_shunt_admittance(branch)
        
        # Handle transformer
        tap = branch.tap_ratio
        if tap == 0.0
            tap = 1.0
        end
        
        shift_rad = branch.phase_shift * π / 180
        tap_complex = tap * exp(im * shift_rad)
        
        # Fill Y-bus matrix
        if tap == 1.0 && branch.phase_shift == 0.0  # Regular line
            Ybus[from_idx, from_idx] += y + b_sh/2
            Ybus[to_idx, to_idx] += y + b_sh/2
            Ybus[from_idx, to_idx] -= y
            Ybus[to_idx, from_idx] -= y
        else  # Transformer
            Ybus[from_idx, from_idx] += y / (tap_complex * conj(tap_complex))
            Ybus[to_idx, to_idx] += y
            Ybus[from_idx, to_idx] -= y / conj(tap_complex)
            Ybus[to_idx, from_idx] -= y / tap_complex
        end
    end
    
    # Add bus shunt admittances
    for (id, bus) in sys.buses
        idx = bus_to_idx[id]
        shunt_y = complex(bus.shunt_g, bus.shunt_b)
        Ybus[idx, idx] += shunt_y
    end
    
    return Ybus, bus_to_idx
end
