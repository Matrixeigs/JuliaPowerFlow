using Polynomials
"""
Comprehensive analysis of simulation results

Parameters:
- sol: ODE solution object

Returns:
- results: dictionary containing all analysis results
"""
function analyze_simulation_results(sol)
    println("  - Extracting time series data...")
    println("  - Calculating performance metrics...")
    println("  - Analyzing stability...")
    println("  - Computing power quality indices...")
    println("  - Evaluating protection performance...")
    
    results = Dict()
    
    # Extract time and state vectors
    t = sol.t
    x = hcat(sol.u...)  # Convert to matrix form
    
    results[:time] = t
    results[:states] = x
    
    # ========================================================================
    # TIME SERIES EXTRACTION
    # ========================================================================
    
    # Generator variables
    results[:gen_id] = x[STATE_INDICES[:gen_id], :]
    results[:gen_iq] = x[STATE_INDICES[:gen_iq], :]
    results[:gen_omega] = x[STATE_INDICES[:gen_ωr], :]
    results[:gen_delta] = x[STATE_INDICES[:gen_δ], :]
    
    # PV system variables
    results[:pv_iL] = x[STATE_INDICES[:pv_iL], :]
    results[:pv_Vdc] = x[STATE_INDICES[:pv_Vdc], :]
    
    # BESS variables
    results[:bess_SOC] = x[STATE_INDICES[:bess_SOC], :]
    results[:bess_Vdc] = x[STATE_INDICES[:bess_Vdc], :]
    results[:bess_ibat] = x[STATE_INDICES[:bess_ibat], :]
    
    # Motor variables
    results[:motor_omega] = x[STATE_INDICES[:motor_ωr], :]
    
    # EV variables
    results[:ev_SOC] = x[STATE_INDICES[:ev_SOC], :]
    results[:ev_iac] = x[STATE_INDICES[:ev_iac], :]
    
    # Network variables
    results[:bus_va] = x[STATE_INDICES[:bus_va], :]
    results[:bus_vb] = x[STATE_INDICES[:bus_vb], :]
    results[:bus_vc] = x[STATE_INDICES[:bus_vc], :]
    results[:line_ia] = x[STATE_INDICES[:line_ia], :]
    results[:line_ib] = x[STATE_INDICES[:line_ib], :]
    results[:line_ic] = x[STATE_INDICES[:line_ic], :]
    results[:fault_current] = x[STATE_INDICES[:fault_current], :]
    
    # ========================================================================
    # CALCULATED QUANTITIES
    # ========================================================================
    
    # System frequency
    results[:frequency] = results[:gen_omega] * F_BASE
    
    # Bus voltage magnitudes
    results[:bus_V_mag] = sqrt.(results[:bus_va].^2 + results[:bus_vb].^2 + results[:bus_vc].^2) / √3
    
    # Line current magnitudes
    results[:line_I_mag] = sqrt.(results[:line_ia].^2 + results[:line_ib].^2 + results[:line_ic].^2) / √3
    
    # Generator power
    results[:gen_P] = results[:gen_id] .* results[:bus_va] + results[:gen_iq] .* results[:bus_vb]
    results[:gen_Q] = results[:gen_iq] .* results[:bus_va] - results[:gen_id] .* results[:bus_vb]
    
    # PV power
    results[:pv_P] = results[:pv_iL] .* results[:pv_Vdc] * 0.95  # Include inverter efficiency
    
    # BESS power
    results[:bess_P] = results[:bess_ibat] .* results[:bess_Vdc] / 1000  # kW
    
    # Motor power
    results[:motor_P] = results[:motor_omega] .* 50.0 / 1000  # Approximate mechanical power in kW
    
    # EV charging power
    results[:ev_P] = results[:ev_iac] .* 240.0 / 1000  # kW
    
    # ========================================================================
    # PERFORMANCE METRICS
    # ========================================================================
    
    metrics = Dict()
    
    # Frequency performance
    freq_deviation = abs.(results[:frequency] .- F_BASE)
    metrics[:max_freq_deviation] = maximum(freq_deviation)
    metrics[:avg_freq_deviation] = mean(freq_deviation)
    metrics[:freq_nadir] = minimum(results[:frequency])
    metrics[:freq_zenith] = maximum(results[:frequency])
    
    # Voltage performance
    voltage_deviation = abs.(results[:bus_V_mag] .- 1.0)
    metrics[:max_voltage_deviation] = maximum(voltage_deviation)
    metrics[:avg_voltage_deviation] = mean(voltage_deviation)
    metrics[:min_voltage] = minimum(results[:bus_V_mag])
    metrics[:max_voltage] = maximum(results[:bus_V_mag])
    
    # Generator performance
    metrics[:max_gen_current] = maximum(sqrt.(results[:gen_id].^2 + results[:gen_iq].^2))
    metrics[:max_power_angle] = maximum(abs.(results[:gen_delta]))
    metrics[:gen_stability_margin] = π/2 - metrics[:max_power_angle]
    
    # BESS performance
    metrics[:bess_energy_throughput] = sum(abs.(results[:bess_P])) * DT / 3600  # kWh
    metrics[:bess_SOC_range] = maximum(results[:bess_SOC]) - minimum(results[:bess_SOC])
    metrics[:bess_efficiency] = calculate_bess_efficiency(results[:bess_P], results[:bess_SOC])
    
    # PV performance
    metrics[:pv_energy_production] = sum(max.(results[:pv_P], 0)) * DT / 3600  # kWh
    metrics[:pv_capacity_factor] = mean(results[:pv_P]) / maximum(results[:pv_P])
    
    # Motor performance
    motor_speed_deviation = abs.(results[:motor_omega] .- mean(results[:motor_omega]))
    metrics[:max_motor_speed_deviation] = maximum(motor_speed_deviation)
    metrics[:motor_efficiency] = calculate_motor_efficiency(results[:motor_P], results[:motor_omega])
    
    # EV charging performance
    metrics[:ev_energy_charged] = sum(max.(results[:ev_P], 0)) * DT / 3600  # kWh
    metrics[:ev_charging_efficiency] = calculate_ev_efficiency(results[:ev_P], results[:ev_SOC])
    
    results[:metrics] = metrics
    
    # ========================================================================
    # STABILITY ANALYSIS
    # ========================================================================
    
    stability = Dict()
    
    # Small-signal stability (eigenvalue analysis approximation)
    stability[:frequency_oscillations] = analyze_oscillations(results[:frequency], t)
    stability[:voltage_oscillations] = analyze_oscillations(results[:bus_V_mag], t)
    stability[:power_oscillations] = analyze_oscillations(results[:gen_P], t)
    
    # Transient stability
    stability[:critical_clearing_time] = estimate_critical_clearing_time(results[:gen_delta], t)
    stability[:first_swing_stability] = check_first_swing_stability(results[:gen_delta])
    
    # Voltage stability
    stability[:voltage_stability_margin] = calculate_voltage_stability_margin(results[:bus_V_mag], results[:gen_Q])
    
    results[:stability] = stability
    
    # ========================================================================
    # POWER QUALITY ANALYSIS
    # ========================================================================
    
    power_quality = Dict()
    
    # Voltage unbalance
    power_quality[:voltage_unbalance] = calculate_voltage_unbalance_time_series(
        results[:bus_va], results[:bus_vb], results[:bus_vc]
    )
    
    # Harmonic analysis (simplified - would need FFT for full analysis)
    power_quality[:voltage_thd] = estimate_thd(results[:bus_va], t)
    power_quality[:current_thd] = estimate_thd(results[:line_ia], t)
    
    # Flicker analysis
    power_quality[:voltage_flicker] = calculate_flicker_index(results[:bus_V_mag], t)
    
    # Frequency quality
    power_quality[:rocof] = calculate_rocof(results[:frequency], t)  # Rate of Change of Frequency
    
    results[:power_quality] = power_quality
    
    # ========================================================================
    # PROTECTION ANALYSIS
    # ========================================================================
    
    protection = Dict()
    
    # Fault detection performance
    fault_start_idx = findfirst(t .>= FAULT_PARAMS[:t_start])
    fault_end_idx = findfirst(t .>= FAULT_PARAMS[:t_end])
    
    if fault_start_idx !== nothing && fault_end_idx !== nothing
        protection[:fault_current_peak] = maximum(results[:fault_current][fault_start_idx:fault_end_idx])
        protection[:fault_duration] = FAULT_PARAMS[:t_end] - FAULT_PARAMS[:t_start]
        protection[:voltage_dip] = 1.0 - minimum(results[:bus_V_mag][fault_start_idx:fault_end_idx])
    end
    
    # Protection relay performance
    protection[:overcurrent_operations] = count_protection_operations(results[:line_I_mag], 1.5)
    protection[:undervoltage_operations] = count_protection_operations(results[:bus_V_mag], 0.8, :under)
    protection[:underfrequency_operations] = count_protection_operations(results[:frequency], 59.0, :under)
    
    results[:protection] = protection
    
    # ========================================================================
    # ECONOMIC ANALYSIS
    # ========================================================================
    
    economics = Dict()
    
    # Energy costs and revenues (simplified)
    electricity_price = 0.10  # $/kWh
    
    economics[:gen_revenue] = sum(max.(results[:gen_P], 0)) * DT / 3600 * electricity_price
    economics[:pv_revenue] = sum(max.(results[:pv_P], 0)) * DT / 3600 * electricity_price
    economics[:bess_cost] = sum(max.(-results[:bess_P], 0)) * DT / 3600 * electricity_price
    economics[:bess_revenue] = sum(max.(results[:bess_P], 0)) * DT / 3600 * electricity_price
    economics[:ev_cost] = sum(max.(results[:ev_P], 0)) * DT / 3600 * electricity_price
    
    economics[:net_cost] = economics[:gen_revenue] + economics[:pv_revenue] + 
                          economics[:bess_revenue] - economics[:bess_cost] - economics[:ev_cost]
    
    results[:economics] = economics
    
    # ========================================================================
    # ENVIRONMENTAL ANALYSIS
    # ========================================================================
    
    environmental = Dict()
    
    # CO2 emissions (simplified)
    gen_emission_factor = 0.5  # kg CO2/kWh for natural gas
    grid_emission_factor = 0.4  # kg CO2/kWh average grid
    
    environmental[:gen_emissions] = sum(max.(results[:gen_P], 0)) * DT / 3600 * gen_emission_factor
    environmental[:pv_emissions_avoided] = sum(max.(results[:pv_P], 0)) * DT / 3600 * grid_emission_factor
    environmental[:net_emissions] = environmental[:gen_emissions] - environmental[:pv_emissions_avoided]
    
    results[:environmental] = environmental
    
    println("  - Analysis completed successfully!")
    
    return results
end

"""
Calculate BESS efficiency based on power and SOC data
"""
function calculate_bess_efficiency(P_bess::Vector{Float64}, SOC::Vector{Float64})
    # Separate charging and discharging periods
    charging_mask = P_bess .< 0
    discharging_mask = P_bess .> 0
    
    if sum(charging_mask) > 0 && sum(discharging_mask) > 0
        energy_in = sum(abs.(P_bess[charging_mask])) * DT / 3600
        energy_out = sum(P_bess[discharging_mask]) * DT / 3600
        
        if energy_in > 0
            return energy_out / energy_in
        end
    end
    
    return 0.9  # Default efficiency
end

"""
Calculate motor efficiency
"""
function calculate_motor_efficiency(P_motor::Vector{Float64}, omega::Vector{Float64})
    # Motor efficiency varies with load and speed
    rated_speed = ω_BASE
    rated_power = 50.0  # kW
    
    # Normalized speed and power
    speed_ratio = mean(omega) / rated_speed
    power_ratio = mean(abs.(P_motor)) / rated_power
    
    # Simplified efficiency curve
    if power_ratio > 0.1
        efficiency = 0.85 + 0.1 * power_ratio * (1 - power_ratio^2) * speed_ratio
        return min(0.95, max(0.7, efficiency))
    else
        return 0.7  # Low efficiency at very low loads
    end
end

"""
Calculate EV charging efficiency
"""
function calculate_ev_efficiency(P_ev::Vector{Float64}, SOC::Vector{Float64})
    # EV charging efficiency typically decreases as SOC increases
    avg_soc = mean(SOC)
    
    if avg_soc < 0.8
        return 0.92  # High efficiency for normal charging
    else
        return 0.85  # Lower efficiency for high SOC charging
    end
end

"""
Analyze oscillations in a signal
"""
function analyze_oscillations(signal::Vector{Float64}, t::Vector{Float64})
    oscillation_data = Dict()
    
    # Remove DC component
    signal_ac = signal .- mean(signal)
    
    # Calculate RMS of oscillations
    oscillation_data[:rms] = sqrt(mean(signal_ac.^2))
    
    # Peak-to-peak amplitude
    oscillation_data[:peak_to_peak] = maximum(signal) - minimum(signal)
    
    # Estimate dominant frequency (simplified)
    if length(signal) > 100
        # Find zero crossings
        zero_crossings = []
        for i in 2:length(signal_ac)
            if signal_ac[i-1] * signal_ac[i] < 0
                push!(zero_crossings, i)
            end
        end
        
        if length(zero_crossings) > 2
            # Estimate frequency from zero crossings
            avg_period = 2 * mean(diff(zero_crossings)) * (t[2] - t[1])
            oscillation_data[:dominant_frequency] = 1.0 / avg_period
        else
            oscillation_data[:dominant_frequency] = 0.0
        end
    else
        oscillation_data[:dominant_frequency] = 0.0
    end
    
    # Damping estimation (exponential decay fit)
    if oscillation_data[:rms] > 1e-6
        envelope = abs.(signal_ac)
        if length(envelope) > 50
            # Fit exponential decay to envelope
            t_fit = t[1:min(length(t), 100)]
            env_fit = envelope[1:length(t_fit)]
            
            # Simple linear fit to log(envelope)
            log_env = log.(max.(env_fit, 1e-10))
            if !any(isinf.(log_env))
                p = polyfit(t_fit, log_env, 1)
                oscillation_data[:damping_ratio] = -p[1] / (2 * π * oscillation_data[:dominant_frequency])
            else
                oscillation_data[:damping_ratio] = 0.1
            end
        else
            oscillation_data[:damping_ratio] = 0.1
        end
    else
        oscillation_data[:damping_ratio] = 1.0  # Critically damped
    end
    
    return oscillation_data
end

"""
Estimate critical clearing time for transient stability
"""
function estimate_critical_clearing_time(delta::Vector{Float64}, t::Vector{Float64})
    # Find maximum power angle during fault
    fault_start_idx = findfirst(t .>= FAULT_PARAMS[:t_start])
    fault_end_idx = findfirst(t .>= FAULT_PARAMS[:t_end])
    
    if fault_start_idx !== nothing && fault_end_idx !== nothing
        max_delta_during_fault = maximum(abs.(delta[fault_start_idx:fault_end_idx]))
        
        # Critical angle is typically around 120-150 degrees
        critical_angle = 2.0  # radians (≈ 115 degrees)
        
        if max_delta_during_fault < critical_angle
            # System remained stable, CCT is longer than actual clearing time
            return FAULT_PARAMS[:t_end] - FAULT_PARAMS[:t_start] + 0.1
        else
            # System became unstable, CCT is shorter than actual clearing time
            return max(0.05, FAULT_PARAMS[:t_end] - FAULT_PARAMS[:t_start] - 0.05)
        end
    else
        return 0.2  # Default CCT estimate
    end
end

"""
Check first swing stability
"""
function check_first_swing_stability(delta::Vector{Float64})
    max_delta = maximum(abs.(delta))
    
    # First swing stability criteria
    if max_delta < π/2  # 90 degrees
        return :stable
    elseif max_delta < 2π/3  # 120 degrees
        return :marginal
    else
        return :unstable
    end
end

"""
Calculate voltage stability margin
"""
function calculate_voltage_stability_margin(V_mag::Vector{Float64}, Q_gen::Vector{Float64})
    # Simplified P-V curve analysis
    min_voltage = minimum(V_mag)
    
    # Voltage stability margin based on minimum voltage
    if min_voltage > 0.95
        return :high_margin
    elseif min_voltage > 0.90
        return :adequate_margin
    elseif min_voltage > 0.85
        return :low_margin
    else
        return :critical_margin
    end
end

"""
Calculate voltage unbalance time series
"""
function calculate_voltage_unbalance_time_series(Va::Vector{Float64}, Vb::Vector{Float64}, Vc::Vector{Float64})
    unbalance = zeros(length(Va))
    
    for i in 1:length(Va)
        # Convert to complex phasors
        Va_complex = Va[i] + 0im
        Vb_complex = Vb[i] * exp(-im * 2π/3)
        Vc_complex = Vc[i] * exp(im * 2π/3)
        
        # Calculate sequence components
        V0, V1, V2 = symmetrical_components(Va_complex, Vb_complex, Vc_complex)
        
        # Voltage unbalance factor
        if abs(V1) > 1e-6
            unbalance[i] = abs(V2) / abs(V1) * 100.0
        else
            unbalance[i] = 0.0
        end
    end
    
    return unbalance
end

"""
Estimate THD from time series (simplified)
"""
function estimate_thd(signal::Vector{Float64}, t::Vector{Float64})
    # This is a simplified THD estimation
    # Full implementation would require FFT analysis
    
    # Remove DC component
    signal_ac = signal .- mean(signal)
    
    # Calculate fundamental frequency component
    ω_fund = ω_BASE
    dt = t[2] - t[1]
    
    # Correlation with fundamental frequency
    fund_cos = cos.(ω_fund .* t)
    fund_sin = sin.(ω_fund .* t)
    
    a1 = 2 * mean(signal_ac .* fund_cos)
    b1 = 2 * mean(signal_ac .* fund_sin)
    
    fund_magnitude = sqrt(a1^2 + b1^2)
    
    # Total RMS
    total_rms = sqrt(mean(signal_ac.^2))
    
    # THD estimation
    if fund_magnitude > 1e-6
        harmonic_rms = sqrt(max(0, total_rms^2 - (fund_magnitude/√2)^2))
        thd = harmonic_rms / (fund_magnitude/√2) * 100.0
        return min(thd, 50.0)  # Cap at 50%
    else
        return 0.0
    end
end

"""
Calculate flicker index
"""
function calculate_flicker_index(V_mag::Vector{Float64}, t::Vector{Float64})
    # Voltage flicker based on IEC 61000-4-15
    
    # Calculate voltage variations
    V_ref = mean(V_mag)
    dV = abs.(diff(V_mag))
    dt = mean(diff(t))
    
    # Flicker severity (simplified)
    if length(dV) > 0
        flicker_rate = mean(dV) / dt
        flicker_amplitude = std(V_mag) / V_ref
        
        # Simplified flicker index
        flicker_index = flicker_amplitude * sqrt(flicker_rate) * 100
        return min(flicker_index, 10.0)  # Cap at reasonable value
    else
        return 0.0
    end
end

"""
Calculate Rate of Change of Frequency (ROCOF)
"""
function calculate_rocof(frequency::Vector{Float64}, t::Vector{Float64})
    if length(frequency) < 2
        return zeros(length(frequency))
    end
    
    dt = mean(diff(t))
    rocof = diff(frequency) / dt
    
    # Add zero at the beginning to match array size
    return [0.0; rocof]
end

"""
Count protection operations
"""
function count_protection_operations(signal::Vector{Float64}, threshold::Float64, direction::Symbol=:over)
    operations = 0
    in_operation = false
    
    for value in signal
        if direction == :over
            condition = value > threshold
        else  # :under
            condition = value < threshold
        end
        
        if condition && !in_operation
            operations += 1
            in_operation = true
        elseif !condition
            in_operation = false
        end
    end
    
    return operations
end

"""
Generate comprehensive simulation report
"""
function generate_simulation_report(results::Dict, components, filename::String="simulation_report.txt")
    open(filename, "w") do file
        println(file, "="^80)
        println(file, "ELECTROMAGNETIC TRANSIENT SIMULATION REPORT")
        println(file, "="^80)
        println(file, "Generated: $(Dates.now())")
        println(file, "Simulation Duration: $(results[:time][end] - results[:time][1]) seconds")
        println(file, "Time Steps: $(length(results[:time]))")
        println(file, "")
        
        # System Configuration
        println(file, "SYSTEM CONFIGURATION")
        println(file, "-"^40)
        println(file, "Generator Rating: $(components.generator.S_rated) MVA")
        println(file, "PV System Rating: $(components.pv_system.P_rated) kW")
        println(file, "BESS Capacity: $(components.bess.Q_max) Ah")
        println(file, "Motor Rating: $(components.motor.P_rated) kW")
        println(file, "")
        
        # Performance Metrics
        println(file, "PERFORMANCE METRICS")
        println(file, "-"^40)
        metrics = results[:metrics]
        println(file, "Maximum Frequency Deviation: $(round(metrics[:max_freq_deviation], digits=3)) Hz")
        println(file, "Maximum Voltage Deviation: $(round(metrics[:max_voltage_deviation], digits=3)) pu")
        println(file, "Generator Stability Margin: $(round(metrics[:gen_stability_margin], digits=3)) rad")
        println(file, "BESS Energy Throughput: $(round(metrics[:bess_energy_throughput], digits=2)) kWh")
        println(file, "PV Energy Production: $(round(metrics[:pv_energy_production], digits=2)) kWh")
        println(file, "")
        
        # Stability Analysis
        println(file, "STABILITY ANALYSIS")
        println(file, "-"^40)
        stability = results[:stability]
        println(file, "First Swing Stability: $(stability[:first_swing_stability])")
        println(file, "Voltage Stability Margin: $(stability[:voltage_stability_margin])")
        println(file, "Critical Clearing Time: $(round(stability[:critical_clearing_time], digits=3)) s")
        println(file, "")
        
        # Power Quality
        println(file, "POWER QUALITY")
        println(file, "-"^40)
        pq = results[:power_quality]
        println(file, "Voltage THD: $(round(pq[:voltage_thd], digits=2))%")
        println(file, "Current THD: $(round(pq[:current_thd], digits=2))%")
        println(file, "Maximum Voltage Unbalance: $(round(maximum(pq[:voltage_unbalance]), digits=2))%")
        println(file, "Maximum ROCOF: $(round(maximum(abs.(pq[:rocof])), digits=2)) Hz/s")
        println(file, "")
        
        # Protection Performance
        if haskey(results, :protection)
            println(file, "PROTECTION PERFORMANCE")
            println(file, "-"^40)
            protection = results[:protection]
            if haskey(protection, :fault_current_peak)
                println(file, "Fault Current Peak: $(round(protection[:fault_current_peak], digits=2)) pu")
                println(file, "Fault Duration: $(protection[:fault_duration]) s")
                println(file, "Voltage Dip: $(round(protection[:voltage_dip], digits=3)) pu")
            end
            println(file, "Overcurrent Operations: $(protection[:overcurrent_operations])")
            println(file, "Undervoltage Operations: $(protection[:undervoltage_operations])")
            println(file, "")
        end
        
        # Economic Analysis
        println(file, "ECONOMIC ANALYSIS")
        println(file, "-"^40)
        economics = results[:economics]
        println(file, "Generator Revenue: \$$(round(economics[:gen_revenue], digits=2))")
        println(file, "PV Revenue: \$$(round(economics[:pv_revenue], digits=2))")
        println(file, "BESS Net Revenue: \$$(round(economics[:bess_revenue] - economics[:bess_cost], digits=2))")
        println(file, "EV Charging Cost: \$$(round(economics[:ev_cost], digits=2))")
        println(file, "Net System Cost: \$$(round(economics[:net_cost], digits=2))")
        println(file, "")
        
        # Environmental Impact
        println(file, "ENVIRONMENTAL IMPACT")
        println(file, "-"^40)
        env = results[:environmental]
        println(file, "Generator Emissions: $(round(env[:gen_emissions], digits=2)) kg CO2")
        println(file, "PV Emissions Avoided: $(round(env[:pv_emissions_avoided], digits=2)) kg CO2")
        println(file, "Net Emissions: $(round(env[:net_emissions], digits=2)) kg CO2")
        println(file, "")
        
        # Recommendations
        println(file, "RECOMMENDATIONS")
        println(file, "-"^40)
        
        if metrics[:max_freq_deviation] > 0.5
            println(file, "- Consider additional frequency regulation resources")
        end
        
        if metrics[:max_voltage_deviation] > 0.05
            println(file, "- Voltage regulation may need improvement")
        end
        
        if stability[:first_swing_stability] == :marginal
            println(file, "- Transient stability margins are low")
        end
        
        if pq[:voltage_thd] > 5.0
            println(file, "- Harmonic filtering may be required")
        end
        
        if maximum(pq[:voltage_unbalance]) > 2.0
            println(file, "- Address voltage unbalance issues")
        end
        
        println(file, "")
        println(file, "="^80)
    end
    
    println("  - Simulation report generated: $filename")
end

"""
Create visualization plots of simulation results
"""
function create_visualization_plots(results::Dict)
    println("  - Creating visualization plots...")
    
    # This would typically use a plotting package like Plots.jl or PyPlot.jl
    # For now, we'll create the data structures that would be used for plotting
    
    plots_data = Dict()
    
    # Time vector
    t = results[:time]
    
    # Plot 1: System frequency
    plots_data[:frequency_plot] = Dict(
        :x => t,
        :y => results[:frequency],
        :title => "System Frequency",
        :xlabel => "Time (s)",
        :ylabel => "Frequency (Hz)",
        :ylim => [59.5, 60.5]
    )
    
    # Plot 2: Bus voltages
    plots_data[:voltage_plot] = Dict(
        :x => t,
        :y => [results[:bus_va], results[:bus_vb], results[:bus_vc]],
        :labels => ["Phase A", "Phase B", "Phase C"],
        :title => "Bus Voltages",
        :xlabel => "Time (s)",
        :ylabel => "Voltage (pu)"
    )
    
    # Plot 3: Generator power angle
    plots_data[:power_angle_plot] = Dict(
        :x => t,
        :y => results[:gen_delta] * 180/π,  # Convert to degrees
        :title => "Generator Power Angle",
        :xlabel => "Time (s)",
        :ylabel => "Power Angle (degrees)"
    )
    
    # Plot 4: BESS SOC and Power
    plots_data[:bess_plot] = Dict(
        :x => t,
        :y1 => results[:bess_SOC] * 100,  # Percentage
        :y2 => results[:bess_P],
        :title => "BESS Performance",
        :xlabel => "Time (s)",
        :ylabel1 => "SOC (%)",
        :ylabel2 => "Power (kW)"
    )
    
    # Plot 5: Power flows
    plots_data[:power_plot] = Dict(
        :x => t,
        :y => [results[:gen_P], results[:pv_P], results[:bess_P], results[:motor_P]],
        :labels => ["Generator", "PV", "BESS", "Motor"],
        :title => "System Power Flows",
        :xlabel => "Time (s)",
        :ylabel => "Power (kW)"
    )
    
    # Plot 6: Fault current (if fault occurred)
    if maximum(abs.(results[:fault_current])) > 0.1
        plots_data[:fault_current_plot] = Dict(
            :x => t,
            :y => results[:fault_current],
            :title => "Fault Current",
            :xlabel => "Time (s)",
            :ylabel => "Current (pu)"
        )
    end
    
    return plots_data
end