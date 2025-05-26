"""
Generator D-Shape Capability Curve Visualization (Fixed Version)

This script demonstrates the characteristic D-shaped capability region of synchronous generators,
showing how various physical constraints combine to create the operating envelope.
"""

using Plots, LinearAlgebra, Printf

"""
    calculate_underexcitation_q(p, eq_min, ut, xd)

Calculate Q for underexcitation limit at given P.
"""
function calculate_underexcitation_q(p, eq_min, ut, xd)
    if p ≤ 0
        return (eq_min * ut - ut^2) / xd
    end
    
    sin_delta = p * xd / (eq_min * ut)
    if abs(sin_delta) > 1
        return NaN  # Changed from -Inf to NaN for better handling
    end
    
    cos_delta = sqrt(1 - sin_delta^2)
    return (eq_min * ut * cos_delta - ut^2) / xd
end

"""
    calculate_overexcitation_q(p, eq_max, ut, xd)

Calculate Q for overexcitation limit at given P.
"""
function calculate_overexcitation_q(p, eq_max, ut, xd)
    if p ≤ 0
        return (eq_max * ut - ut^2) / xd
    end
    
    sin_delta = p * xd / (eq_max * ut)
    if abs(sin_delta) > 1
        return NaN  # Changed from Inf to NaN for better handling
    end
    
    cos_delta = sqrt(1 - sin_delta^2)
    return (eq_max * ut * cos_delta - ut^2) / xd
end

"""
    construct_d_shape_boundary(p_min, p_max, s_max, eq_min, eq_max, ut, xd, δ_max_rad, q_min)

Construct the complete boundary of the D-shaped capability region.
"""
function construct_d_shape_boundary(p_min, p_max, s_max, eq_min, eq_max, ut, xd, δ_max_rad, q_min)
    boundary_p = Float64[]
    boundary_q = Float64[]
    
    # Start from P_min at Q_min
    push!(boundary_p, p_min)
    push!(boundary_q, q_min)
    
    # Follow Q_min line until it intersects with underexcitation limit
    p_current = p_min
    while p_current ≤ p_max
        # Check if we're still on Q_min or need to follow underexcitation
        q_under = calculate_underexcitation_q(p_current, eq_min, ut, xd)
        
        if !isnan(q_under) && q_under > q_min
            # Switch to underexcitation limit
            break
        end
        
        p_current += 0.01
        push!(boundary_p, p_current)
        push!(boundary_q, q_min)
    end
    
    # Follow underexcitation limit
    p_range = range(p_current, p_max, length=50)
    for p in p_range
        q = calculate_underexcitation_q(p, eq_min, ut, xd)
        if !isnan(q)  # Only add valid points
            push!(boundary_p, p)
            push!(boundary_q, q)
        end
    end
    
    # Follow right boundary (P_max or stability limit)
    p_stab_max = (eq_max * ut / xd) * sin(δ_max_rad)
    p_right = min(p_max, p_stab_max)
    
    if p_right < p_max
        # Follow stability limit
        eq_range = range(eq_max, eq_min, length=50)
        for eq in eq_range
            p = (eq * ut / xd) * sin(δ_max_rad)
            q = (eq * ut * cos(δ_max_rad) - ut^2) / xd
            if p ≥ p_min
                push!(boundary_p, p)
                push!(boundary_q, q)
            end
        end
    else
        # Vertical line at P_max
        q_bottom = calculate_underexcitation_q(p_max, eq_min, ut, xd)
        q_top = calculate_overexcitation_q(p_max, eq_max, ut, xd)
        
        # Handle NaN values
        if isnan(q_bottom)
            q_bottom = q_min
        end
        if isnan(q_top)
            q_top = sqrt(max(0, s_max^2 - p_max^2))
        end
        
        if q_bottom < q_top  # Only create range if valid
            q_range = range(q_bottom, q_top, length=20)
            for q in q_range
                push!(boundary_p, p_max)
                push!(boundary_q, q)
            end
        end
    end
    
    # Follow overexcitation limit back to start
    p_range_rev = range(p_right, p_min, length=50)
    for p in p_range_rev
        q = calculate_overexcitation_q(p, eq_max, ut, xd)
        if !isnan(q)  # Only add valid points
            # Also check thermal limit
            if p^2 + q^2 > s_max^2
                q_thermal = sqrt(max(0, s_max^2 - p^2))
                q = min(q, q_thermal)
            end
            push!(boundary_p, p)
            push!(boundary_q, q)
        end
    end
    
    # Close the boundary
    push!(boundary_p, p_min)
    push!(boundary_q, q_min)
    
    return boundary_p, boundary_q
end

"""
    plot_d_shape_generator(; kwargs...)

Plot the complete D-shaped capability curve for a synchronous generator.
"""
function plot_d_shape_generator(;
    p_max=0.9,           # Maximum active power (p.u.)
    p_min=0.2,           # Minimum active power (p.u.)
    s_max=1.0,           # Maximum apparent power (p.u.)
    eq_min=0.8,          # Minimum internal voltage (p.u.)
    eq_max=1.4,          # Maximum internal voltage (p.u.)
    ut=1.0,              # Terminal voltage (p.u.)
    xd=0.8,              # Direct-axis reactance (p.u.)
    delta_max=60.0,      # Maximum power angle (degrees)
    q_min_constant=-0.3, # Constant Q minimum (p.u.)
    title="Generator D-Shape Capability Curve",
    n_points=200)
    
    # Convert angle to radians
    δ_max_rad = delta_max * π / 180
    
    # Initialize arrays for boundary curves
    P_thermal = Float64[]
    Q_thermal = Float64[]
    P_overexcited = Float64[]
    Q_overexcited = Float64[]
    P_underexcited = Float64[]
    Q_underexcited = Float64[]
    P_stability = Float64[]
    Q_stability = Float64[]
    
    # 1. Thermal (Apparent Power) Limit: P² + Q² = S_max²
    θ = range(-π/2, π/2, length=n_points)
    for angle in θ
        p = s_max * cos(angle)
        q = s_max * sin(angle)
        if p_min ≤ p ≤ p_max  # Only include feasible active power range
            push!(P_thermal, p)
            push!(Q_thermal, q)
        end
    end
    
    # 2. Overexcitation Limit (Field Current Limit)
    p_range = range(p_min, min(p_max, eq_max * ut / xd), length=n_points)
    for p in p_range
        q = calculate_overexcitation_q(p, eq_max, ut, xd)
        if !isnan(q)  # Changed from !isinf(q)
            push!(P_overexcited, p)
            push!(Q_overexcited, q)
        end
    end
    
    # 3. Underexcitation Limit
    for p in p_range
        q = calculate_underexcitation_q(p, eq_min, ut, xd)
        if !isnan(q)  # Changed from !isinf(q)
            push!(P_underexcited, p)
            push!(Q_underexcited, q)
        end
    end
    
    # 4. Stability Limit (Power Angle Limit)
    eq_range = range(eq_min, eq_max, length=n_points)
    for eq in eq_range
        p_stab = (eq * ut / xd) * sin(δ_max_rad)
        if p_min ≤ p_stab ≤ p_max
            q_stab = (eq * ut * cos(δ_max_rad) - ut^2) / xd
            push!(P_stability, p_stab)
            push!(Q_stability, q_stab)
        end
    end
    
    # Create the main plot with better legend positioning
    plot_obj = plot(xlabel="Reactive Power Q (p.u.)", 
                   ylabel="Active Power P (p.u.)",
                   title=title,
                   grid=true,
                   legend=:outertopright,  # Changed legend position to avoid overlap
                   size=(800, 600),       # Increased width to accommodate legend
                   left_margin=5Plots.mm,
                   right_margin=15Plots.mm,  # Extra margin for external legend
                   aspect_ratio=:equal)
    
    # Plot thermal limit (circular arc) - swapped Q and P
    if !isempty(P_thermal)
        plot!(Q_thermal, P_thermal,
              label="Thermal Limit",
              linewidth=2,
              color=:red,
              linestyle=:solid)
    end
    
    # Plot overexcitation limit (upper D boundary) - swapped Q and P
    if !isempty(P_overexcited)
        plot!(Q_overexcited, P_overexcited,
              label="Overexcitation",
              linewidth=2,
              color=:blue,
              linestyle=:dash)
    end
    
    # Plot underexcitation limit (lower D boundary) - swapped Q and P
    if !isempty(P_underexcited)
        plot!(Q_underexcited, P_underexcited,
              label="Underexcitation",
              linewidth=2,
              color=:green,
              linestyle=:dot)
    end
    
    # Plot stability limit (right D boundary) - swapped Q and P
    if !isempty(P_stability)
        plot!(Q_stability, P_stability,
              label="Stability Limit",
              linewidth=2,
              color=:orange,
              linestyle=:dashdot)
    end
    
    # Add active power limits as horizontal lines (now on y-axis)
    hline!([p_min], label="P_min", color=:black, linewidth=1, linestyle=:solid)
    hline!([p_max], label="P_max", color=:black, linewidth=1, linestyle=:solid)
    
    # Add constant Q_min line as vertical line (now on x-axis)
    vline!([q_min_constant], label="Q_min", 
           color=:purple, linewidth=1, linestyle=:solid)
    
    # Fill the feasible region to show the D-shape
    boundary_p, boundary_q = construct_d_shape_boundary(
        p_min, p_max, s_max, eq_min, eq_max, ut, xd, δ_max_rad, q_min_constant)
    
    if !isempty(boundary_p)
        plot!(boundary_q, boundary_p,  # Swapped Q and P for axes
              fillrange=0,
              fillalpha=0.15,
              fillcolor=:lightblue,
              label="Feasible Region",
              linecolor=:black,
              linewidth=1)
    end
    
    # Mark key operating points - swapped coordinates
    scatter!([0.0], [0.0], 
            markersize=6, 
            color=:red, 
            markershape=:circle,
            label="No Load")
    
    if p_max > 0
        scatter!([0.0], [p_max], 
                markersize=6, 
                color=:blue, 
                markershape=:square,
                label="Rated Power")
    end
    
    return plot_obj
end

# Example usage:
plot_d_shape_generator()
# 
# # Or with custom parameters:
# plot_d_shape_generator(
#     p_max=1.2,
#     eq_min=0.7,
#     eq_max=1.5,
#     delta_max=70.0,
#     title="Custom Generator Capability Curve"
# )
