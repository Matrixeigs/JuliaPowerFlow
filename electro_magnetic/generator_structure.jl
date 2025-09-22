using LinearAlgebra, Plots, FFTW

"""
Synchronous Generator Electromagnetic Analysis
- Armature and field structure modeling
- Multiple pole pairs analysis
- MMF (Magnetomotive Force) waveforms
- Direct and quadrature axes transformations
"""

struct SynchronousGeneratorParameters
    # Physical parameters
    stator_slots::Int           # Number of stator slots
    rotor_poles::Int            # Number of rotor poles (must be even)
    pole_pairs::Int             # Number of pole pairs (p = poles/2)
    stator_radius::Float64      # Stator inner radius (m)
    rotor_radius::Float64       # Rotor outer radius (m)
    air_gap::Float64            # Air gap length (m)
    axial_length::Float64       # Machine axial length (m)
    
    # Winding parameters
    turns_per_coil::Int         # Number of turns per coil
    coils_per_phase::Int        # Number of coils per phase
    winding_factor::Float64     # Winding distribution factor
    
    # Electrical parameters
    rated_voltage::Float64      # Rated line-to-line voltage (V)
    rated_current::Float64      # Rated current (A)
    rated_power::Float64        # Rated power (VA)
    rated_frequency::Float64    # Rated frequency (Hz)
    
    # Machine constants
    Ld::Float64                 # Direct axis synchronous inductance (H)
    Lq::Float64                 # Quadrature axis synchronous inductance (H)
    Lf::Float64                 # Field winding inductance (H)
    Maf::Float64                # Mutual inductance between armature and field (H)
    Ra::Float64                 # Armature resistance (Ω)
    Rf::Float64                 # Field resistance (Ω)
end

struct ArmatureStructure
    slot_positions::Vector{Float64}     # Angular positions of slots (rad)
    coil_sides::Vector{Tuple{Int,Int}}  # Coil side connections (slot pairs)
    phase_distribution::Vector{Int}      # Phase assignment for each slot (1,2,3)
    winding_pattern::Matrix{Int}        # Winding pattern matrix [phase x slot]
    slot_pitch::Float64                 # Slot pitch angle (rad)
    coil_pitch::Float64                 # Coil pitch angle (rad)
end

struct FieldStructure
    pole_positions::Vector{Float64}     # Angular positions of poles (rad)
    pole_width::Float64                 # Angular width of each pole (rad)
    field_distribution::Vector{Float64} # Field strength distribution
    excitation_current::Float64         # Field excitation current (A)
end

struct MMFWaveform
    spatial_positions::Vector{Float64}  # Spatial positions (rad)
    mmf_values::Vector{Float64}         # MMF values at each position
    fundamental_amplitude::Float64       # Fundamental component amplitude
    harmonic_content::Vector{Float64}   # Harmonic amplitudes
    harmonic_orders::Vector{Int}        # Harmonic order numbers
end

struct DQAxisSystem
    d_axis_angle::Float64               # Direct axis position (rad)
    q_axis_angle::Float64               # Quadrature axis position (rad)  
    transformation_matrix::Matrix{Float64} # Park transformation matrix
    inverse_matrix::Matrix{Float64}     # Inverse Park transformation
end

function create_synchronous_generator(poles::Int, slots::Int, rated_power::Float64, 
                                    rated_voltage::Float64, rated_frequency::Float64)
    """Create a synchronous generator with specified parameters"""
    
    pole_pairs = poles ÷ 2
    
    # Calculate machine dimensions (typical proportions)
    rated_current = rated_power / (sqrt(3) * rated_voltage)
    stator_radius = sqrt(rated_power / 1e6) * 0.5  # Empirical scaling
    air_gap = stator_radius * 0.005  # Typical 0.5% of radius
    
    # Winding parameters
    turns_per_coil = Int(ceil(rated_voltage / (4.44 * rated_frequency * 0.1)))
    coils_per_phase = slots ÷ 6  # Two slots per coil, three phases
    
    # Machine constants (typical values scaled by power)
    power_pu = rated_power / 1e6  # Per unit based on 1MVA base
    Ld = 1.8 * power_pu / (rated_current^2)
    Lq = 1.2 * power_pu / (rated_current^2)
    Lf = 2.0 * power_pu / (rated_current^2)
    Maf = 1.5 * power_pu / (rated_current^2)
    Ra = 0.01 * rated_voltage / rated_current
    Rf = 50.0 * Ra
    
    return SynchronousGeneratorParameters(
        slots, poles, pole_pairs, stator_radius, stator_radius * 0.95, 
        air_gap, stator_radius * 1.5, turns_per_coil, coils_per_phase, 0.95,
        rated_voltage, rated_current, rated_power, rated_frequency,
        Ld, Lq, Lf, Maf, Ra, Rf
    )
end

function design_armature_structure(gen::SynchronousGeneratorParameters)
    """Design the armature (stator) structure with three-phase windings"""
    
    # Calculate slot positions
    slot_pitch = 2π / gen.stator_slots
    slot_positions = [(i-1) * slot_pitch for i in 1:gen.stator_slots]
    
    # Design three-phase distributed winding
    slots_per_pole = gen.stator_slots ÷ gen.rotor_poles
    slots_per_pole_per_phase = slots_per_pole ÷ 3
    
    # Phase distribution
    phase_distribution = zeros(Int, gen.stator_slots)
    winding_pattern = zeros(Int, 3, gen.stator_slots)
    
    for slot in 1:gen.stator_slots
        # Determine which pole and phase this slot belongs to
        pole_number = ((slot - 1) ÷ slots_per_pole) + 1
        slot_in_pole = ((slot - 1) % slots_per_pole) + 1
        
        # Phase assignment within each pole
        if slot_in_pole <= slots_per_pole_per_phase
            phase = 1  # Phase A
        elseif slot_in_pole <= 2 * slots_per_pole_per_phase
            phase = 2  # Phase B
        else
            phase = 3  # Phase C
        end
        
        # Account for alternating pole polarity
        polarity = (pole_number % 2 == 1) ? 1 : -1
        
        phase_distribution[slot] = phase
        winding_pattern[phase, slot] = polarity
    end
    
    # Generate coil connections
    coil_sides = Tuple{Int,Int}[]
    coil_pitch_slots = Int(round(gen.stator_slots / gen.rotor_poles * 0.8))  # Short pitch
    
    for slot in 1:gen.stator_slots÷2
        coil_end = slot + coil_pitch_slots
        if coil_end <= gen.stator_slots
            push!(coil_sides, (slot, coil_end))
        end
    end
    
    coil_pitch = coil_pitch_slots * slot_pitch
    
    return ArmatureStructure(slot_positions, coil_sides, phase_distribution, 
                           winding_pattern, slot_pitch, coil_pitch)
end

function design_field_structure(gen::SynchronousGeneratorParameters, excitation_current::Float64)
    """Design the field (rotor) structure with salient poles"""
    
    pole_pitch = 2π / gen.rotor_poles
    pole_positions = [(i-1) * pole_pitch for i in 1:gen.rotor_poles]
    
    # Pole width (typically 60-70% of pole pitch)
    pole_width = 0.65 * pole_pitch
    
    # Field distribution - rectangular approximation
    angular_resolution = 360
    angles = range(0, 2π, length=angular_resolution)
    field_distribution = zeros(angular_resolution)
    
    for (i, angle) in enumerate(angles)
        for pole_pos in pole_positions
            # Check if angle is within pole width
            angle_diff = abs(angle - pole_pos)
            angle_diff = min(angle_diff, 2π - angle_diff)  # Account for wraparound
            
            if angle_diff <= pole_width/2
                # Alternating polarity
                pole_index = Int(round(pole_pos / pole_pitch)) + 1
                polarity = (pole_index % 2 == 1) ? 1 : -1
                field_distribution[i] = polarity * excitation_current
                break
            end
        end
    end
    
    return FieldStructure(pole_positions, pole_width, field_distribution, excitation_current)
end

function calculate_armature_mmf(armature::ArmatureStructure, gen::SynchronousGeneratorParameters,
                               phase_currents::Vector{Float64}, time::Float64=0.0)
    """Calculate armature MMF waveform"""
    
    angular_resolution = 360
    angles = range(0, 2π, length=angular_resolution)
    mmf_total = zeros(angular_resolution)
    
    # Calculate MMF contribution from each phase
    for phase in 1:3
        phase_current = phase_currents[phase]
        
        for (i, angle) in enumerate(angles)
            mmf_contribution = 0.0
            
            # Sum contributions from all slots of this phase
            for slot in 1:gen.stator_slots
                if armature.phase_distribution[slot] == phase
                    slot_angle = armature.slot_positions[slot]
                    winding_polarity = armature.winding_pattern[phase, slot]
                    
                    # MMF contribution from this slot
                    turns = gen.turns_per_coil * gen.winding_factor
                    mmf_at_angle = turns * phase_current * winding_polarity
                    
                    # Spatial distribution (rectangular approximation)
                    angle_diff = abs(angle - slot_angle)
                    angle_diff = min(angle_diff, 2π - angle_diff)
                    
                    if angle_diff <= armature.slot_pitch/2
                        mmf_contribution += mmf_at_angle
                    end
                end
            end
            
            mmf_total[i] += mmf_contribution
        end
    end
    
    # Calculate harmonics using FFT
    mmf_fft = fft(mmf_total)
    fundamental_amplitude = 2 * abs(mmf_fft[gen.pole_pairs + 1]) / angular_resolution
    
    # Extract significant harmonics
    harmonic_orders = Int[]
    harmonic_amplitudes = Float64[]
    
    for h in 1:min(20, angular_resolution÷2)
        amplitude = 2 * abs(mmf_fft[h + 1]) / angular_resolution
        if amplitude > 0.01 * fundamental_amplitude  # Only significant harmonics
            push!(harmonic_orders, h)
            push!(harmonic_amplitudes, amplitude)
        end
    end
    
    return MMFWaveform(collect(angles), mmf_total, fundamental_amplitude, 
                      harmonic_amplitudes, harmonic_orders)
end

function calculate_field_mmf(field::FieldStructure, rotor_angle::Float64=0.0)
    """Calculate field MMF waveform"""
    
    angular_resolution = 360
    angles = range(0, 2π, length=angular_resolution)
    mmf_field = zeros(angular_resolution)
    
    for (i, angle) in enumerate(angles)
        # Rotate field distribution by rotor angle
        rotated_angle = angle - rotor_angle
        rotated_angle = mod(rotated_angle, 2π)
        
        # Interpolate field distribution
        index = Int(round(rotated_angle / (2π) * length(field.field_distribution))) + 1
        index = max(1, min(index, length(field.field_distribution)))
        
        mmf_field[i] = field.field_distribution[index]
    end
    
    # Calculate harmonics
    mmf_fft = fft(mmf_field)
    fundamental_amplitude = 2 * abs(mmf_fft[2]) / angular_resolution
    
    harmonic_orders = Int[]
    harmonic_amplitudes = Float64[]
    
    for h in 1:min(20, angular_resolution÷2)
        amplitude = 2 * abs(mmf_fft[h + 1]) / angular_resolution
        if amplitude > 0.01 * fundamental_amplitude
            push!(harmonic_orders, h)
            push!(harmonic_amplitudes, amplitude)
        end
    end
    
    return MMFWaveform(collect(angles), mmf_field, fundamental_amplitude,
                      harmonic_amplitudes, harmonic_orders)
end

function create_dq_axis_system(rotor_angle::Float64, pole_pairs::Int)
    """Create direct and quadrature axis system"""
    
    # Direct axis aligned with rotor north pole
    d_axis_angle = rotor_angle
    
    # Quadrature axis leads direct axis by 90 electrical degrees
    q_axis_angle = d_axis_angle + π/(2*pole_pairs)
    
    # Park transformation matrix (electrical angle)
    θ_elec = rotor_angle * pole_pairs
    
    # Transformation from abc to dq0
    transformation_matrix = (2/3) * [
        cos(θ_elec)           cos(θ_elec - 2π/3)     cos(θ_elec + 2π/3);
        -sin(θ_elec)         -sin(θ_elec - 2π/3)    -sin(θ_elec + 2π/3);
        1/2                   1/2                    1/2
    ]
    
    # Inverse transformation matrix
    inverse_matrix = [
        cos(θ_elec)          -sin(θ_elec)           1;
        cos(θ_elec - 2π/3)   -sin(θ_elec - 2π/3)   1;
        cos(θ_elec + 2π/3)   -sin(θ_elec + 2π/3)   1
    ]
    
    return DQAxisSystem(d_axis_angle, q_axis_angle, transformation_matrix, inverse_matrix)
end

function abc_to_dq(abc_values::Vector{Float64}, dq_system::DQAxisSystem)
    """Transform three-phase abc quantities to dq0"""
    return dq_system.transformation_matrix * abc_values
end

function dq_to_abc(dq0_values::Vector{Float64}, dq_system::DQAxisSystem)
    """Transform dq0 quantities to three-phase abc"""
    return dq_system.inverse_matrix * dq0_values
end

function plot_generator_structure(gen::SynchronousGeneratorParameters, 
                                armature::ArmatureStructure, field::FieldStructure)
    """Plot the generator cross-sectional structure"""
    
    # Create polar plot of generator structure
    θ = range(0, 2π, length=1000)
    
    # Stator outline
    stator_outer_r = gen.stator_radius * 1.2
    stator_inner_r = gen.stator_radius
    
    # Rotor outline  
    rotor_outer_r = gen.rotor_radius
    
    p1 = plot(title="Synchronous Generator Cross Section", 
              xlabel="Position", ylabel="Position", aspect_ratio=:equal)
    
    # Draw stator
    stator_x_outer = stator_outer_r * cos.(θ)
    stator_y_outer = stator_outer_r * sin.(θ)
    stator_x_inner = stator_inner_r * cos.(θ)
    stator_y_inner = stator_inner_r * sin.(θ)
    
    plot!(p1, stator_x_outer, stator_y_outer, color=:blue, linewidth=2, label="Stator Outer")
    plot!(p1, stator_x_inner, stator_y_inner, color=:blue, linewidth=2, label="Stator Inner")
    
    # Draw rotor
    rotor_x = rotor_outer_r * cos.(θ)
    rotor_y = rotor_outer_r * sin.(θ)
    plot!(p1, rotor_x, rotor_y, color=:red, linewidth=2, label="Rotor")
    
    # Mark stator slots
    for (i, slot_pos) in enumerate(armature.slot_positions)
        x_slot = stator_inner_r * cos(slot_pos)
        y_slot = stator_inner_r * sin(slot_pos)
        phase_color = [:red, :yellow, :blue][armature.phase_distribution[i]]
        scatter!(p1, [x_slot], [y_slot], color=phase_color, markersize=3, 
                label=(i==1 ? "Phase A" : (i==2 ? "Phase B" : (i==3 ? "Phase C" : ""))))
    end
    
    # Mark rotor poles
    for pole_pos in field.pole_positions
        x_pole = rotor_outer_r * 0.8 * cos(pole_pos)
        y_pole = rotor_outer_r * 0.8 * sin(pole_pos)
        scatter!(p1, [x_pole], [y_pole], color=:black, marker=:square, 
                markersize=5, label="Field Poles")
    end
    
    return p1
end

function plot_mmf_waveforms(armature_mmf::MMFWaveform, field_mmf::MMFWaveform, 
                           gen::SynchronousGeneratorParameters)
    """Plot MMF waveforms and their harmonics"""
    
    # Convert angles to electrical degrees
    elec_angles = armature_mmf.spatial_positions * 180/π * gen.pole_pairs
    
    p1 = plot(title="MMF Waveforms", xlabel="Electrical Angle (degrees)", 
              ylabel="MMF (AT)", linewidth=2)
    
    plot!(p1, elec_angles, armature_mmf.mmf_values, label="Armature MMF", color=:blue)
    plot!(p1, elec_angles, field_mmf.mmf_values, label="Field MMF", color=:red)
    plot!(p1, elec_angles, armature_mmf.mmf_values + field_mmf.mmf_values, 
          label="Resultant MMF", color=:green, linestyle=:dash)
    
    # Harmonic content plot
    p2 = plot(title="MMF Harmonic Content", xlabel="Harmonic Order", 
              ylabel="Amplitude (AT)", linewidth=2)
    
    if length(armature_mmf.harmonic_orders) > 0
        bar!(p2, armature_mmf.harmonic_orders, armature_mmf.harmonic_content, 
             alpha=0.7, label="Armature Harmonics", color=:blue)
    end
    
    if length(field_mmf.harmonic_orders) > 0
        bar!(p2, field_mmf.harmonic_orders, field_mmf.harmonic_content, 
             alpha=0.7, label="Field Harmonics", color=:red)
    end
    
    return plot(p1, p2, layout=(2,1), size=(800, 600))
end

function plot_dq_axes(gen::SynchronousGeneratorParameters, dq_system::DQAxisSystem, 
                     rotor_angle::Float64)
    """Plot direct and quadrature axes"""
    
    p = plot(title="Direct and Quadrature Axes", xlabel="Position", ylabel="Position", 
             aspect_ratio=:equal, xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
    
    # Draw stator circle
    θ = range(0, 2π, length=100)
    plot!(p, cos.(θ), sin.(θ), color=:gray, linewidth=2, label="Stator")
    
    # Draw rotor position
    rotor_x = 0.8 * cos(rotor_angle)
    rotor_y = 0.8 * sin(rotor_angle)
    plot!(p, [0, rotor_x], [0, rotor_y], color=:red, linewidth=3, 
          arrow=true, label="Rotor Position")
    
    # Draw d-axis (aligned with rotor)
    d_x = 1.2 * cos(dq_system.d_axis_angle)
    d_y = 1.2 * sin(dq_system.d_axis_angle)
    plot!(p, [0, d_x], [0, d_y], color=:blue, linewidth=2, 
          arrow=true, label="d-axis (Direct)")
    
    # Draw q-axis (90° ahead of d-axis)
    q_x = 1.2 * cos(dq_system.q_axis_angle)
    q_y = 1.2 * sin(dq_system.q_axis_angle)
    plot!(p, [0, q_x], [0, q_y], color=:green, linewidth=2, 
          arrow=true, label="q-axis (Quadrature)")
    
    # Add phase axes for reference
    phase_colors = [:red, :yellow, :blue]
    phase_labels = ["Phase A", "Phase B", "Phase C"]
    for i in 1:3
        phase_angle = (i-1) * 2π/3
        phase_x = cos(phase_angle)
        phase_y = sin(phase_angle)
        plot!(p, [0, phase_x], [0, phase_y], color=phase_colors[i], 
              linewidth=1, linestyle=:dash, label=phase_labels[i])
    end
    
    return p
end

# Example usage and demonstration
function demonstrate_synchronous_generator()
    """Demonstrate synchronous generator electromagnetic analysis"""
    
    println("=== Synchronous Generator Electromagnetic Analysis ===")
    
    # Create a 4-pole, 36-slot synchronous generator
    gen = create_synchronous_generator(4, 36, 1e6, 13800.0, 60.0)  # Fixed: use Float64 for voltage and frequency
    println("Generator: $(gen.rotor_poles) poles, $(gen.stator_slots) slots")
    println("Rated: $(gen.rated_power/1e6) MVA, $(gen.rated_voltage) V, $(gen.rated_frequency) Hz")
    
    # Design armature structure
    armature = design_armature_structure(gen)
    println("Armature designed with $(length(armature.coil_sides)) coils")
    
    # Design field structure
    field_current = 100.0  # A
    field = design_field_structure(gen, field_current)
    println("Field structure with $(length(field.pole_positions)) poles")
    
    # Calculate MMF waveforms
    phase_currents = [1000.0, 1000.0*cos(2π/3), 1000.0*cos(4π/3)]  # Balanced currents
    armature_mmf = calculate_armature_mmf(armature, gen, phase_currents)
    
    rotor_angle = π/4  # 45 degrees
    field_mmf = calculate_field_mmf(field, rotor_angle)
    
    println("Armature MMF fundamental: $(round(armature_mmf.fundamental_amplitude)) AT")
    println("Field MMF fundamental: $(round(field_mmf.fundamental_amplitude)) AT")
    
    # Create dq-axis system
    dq_system = create_dq_axis_system(rotor_angle, gen.pole_pairs)
    
    # Transform currents to dq
    dq_currents = abc_to_dq(phase_currents, dq_system)
    println("dq currents: Id=$(round(dq_currents[1],digits=1)) A, Iq=$(round(dq_currents[2],digits=1)) A")
    
    # Generate plots
    p1 = plot_generator_structure(gen, armature, field)
    p2 = plot_mmf_waveforms(armature_mmf, field_mmf, gen)
    p3 = plot_dq_axes(gen, dq_system, rotor_angle)
    
    # Display plots
    display(plot(p1, p3, layout=(1,2), size=(1200, 500)))
    display(p2)
    
    println("\nAnalysis complete!")
    
    return gen, armature, field, armature_mmf, field_mmf, dq_system
end

println("Synchronous Generator Electromagnetic Analysis loaded!")
println("Run demonstrate_synchronous_generator() for a complete example")
