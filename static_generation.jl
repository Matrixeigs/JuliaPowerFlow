using Plots
using LinearAlgebra

"""
    plot_pq_capability_curve(model::Symbol; kwargs...)

Plot the PQ capability curve of a generator according to the specified model.
Available models: :simple_box, :cylindrical_rotor, :with_delta_limits, :full_constraints

Parameters:
- `model`: Symbol representing the generator model
- `p_max`: Maximum active power output
- `p_min`: Minimum active power output
- `q_max`: Maximum reactive power output
- `q_min`: Minimum reactive power output
- `s_max`: Maximum apparent power
- `delta_max`: Maximum power angle (in degrees)
- `eq_min`: Minimum internal voltage magnitude
- `eq_max`: Maximum internal voltage magnitude
- `additional_args...`: Model-specific parameters
"""
function plot_pq_capability_curve(model::Symbol; kwargs...)
    if model == :simple_box
        return plot_simple_box_model(; kwargs...)
    elseif model == :cylindrical_rotor
        return plot_cylindrical_rotor_model(; kwargs...)
    elseif model == :with_delta_limits
        return plot_with_delta_limits(; kwargs...)
    elseif model == :full_constraints
        return plot_full_constraints(; kwargs...)
    else
        error("Unknown model: $model")
    end
end

"""
    plot_simple_box_model(; p_min=-1.0, p_max=1.0, q_min=-1.0, q_max=1.0, title="Simple Box Model")

Plot the simple box model feasibility region with Q on x-axis and P on y-axis.
"""
function plot_simple_box_model(; p_min=-1.0, p_max=1.0, q_min=-1.0, q_max=1.0, title="Simple Box Model")
    # Create plot with Q on x-axis and P on y-axis
    plt = plot(
        xlims=(q_min*1.1, q_max*1.1), 
        ylims=(p_min*1.1, p_max*1.1),
        xlabel="Q (pu)", 
        ylabel="P (pu)",
        title=title,
        aspect_ratio=:equal,
        legend=false,
        grid=true
    )
    
    # Add y-axis line at x = 0
    plot!(plt, [0, 0], [p_min*1.1, p_max*1.1], linecolor=:black, linewidth=1, alpha=0.5)
    
    # Plot rectangular region (note: Q on x-axis, P on y-axis)
    q_values = [q_min, q_max, q_max, q_min, q_min]
    p_values = [p_min, p_min, p_max, p_max, p_min]
    plot!(plt, q_values, p_values, linecolor=:blue, linewidth=2, fill=true, fillalpha=0.2, fillcolor=:lightblue)
    
    # Add origin point
    scatter!(plt, [0], [0], markersize=4, color=:black)
    
    return plt
end

"""
    plot_cylindrical_rotor_model(; p_min=-1.0, p_max=1.0, s_max=1.0, title="Cylindrical Rotor Model")

Plot the cylindrical rotor model capability curve with apparent power constraint.
Q on x-axis, P on y-axis.
"""
function plot_cylindrical_rotor_model(; p_min=-1.0, p_max=1.0, s_max=1.0, title="Cylindrical Rotor Model")
    # Create plot with Q on x-axis and P on y-axis
    plt = plot(
        xlims=(-s_max*1.1, s_max*1.1),
        ylims=(p_min*1.1, p_max*1.1),
        xlabel="Q (pu)",
        ylabel="P (pu)",
        title=title,
        aspect_ratio=:equal,
        legend=false,
        grid=true
    )
    
    # Add y-axis line at x = 0
    plot!(plt, [0, 0], [p_min*1.1, p_max*1.1], linecolor=:black, linewidth=1, alpha=0.5)
    
    # Plot circular region (S² = P² + Q²)
    θ = range(0, 2π, length=100)
    q_circle = s_max .* cos.(θ)
    p_circle = s_max .* sin.(θ)
    
    # Apply P constraints
    valid_points = findall((p_circle .<= p_max) .& (p_circle .>= p_min))
    
    # Plot the circular region
    plot!(plt, q_circle[valid_points], p_circle[valid_points], linecolor=:blue, linewidth=2)
    
    # Add horizontal lines at P = P_max and P = P_min
    if abs(p_max) < s_max
        q_at_pmax = sqrt(s_max^2 - p_max^2)
        plot!(plt, [-q_at_pmax, q_at_pmax], [p_max, p_max], linecolor=:blue, linewidth=2)
    end
    
    if abs(p_min) < s_max
        q_at_pmin = sqrt(s_max^2 - p_min^2)
        plot!(plt, [-q_at_pmin, q_at_pmin], [p_min, p_min], linecolor=:blue, linewidth=2)
    end
    
    # Fill the region
    filled_points_q = []
    filled_points_p = []
    
    # Add points along P_min
    if abs(p_min) < s_max
        push!(filled_points_q, -q_at_pmin, q_at_pmin)
        push!(filled_points_p, p_min, p_min)
    end
    
    # Add points from the circle
    append!(filled_points_q, q_circle[valid_points])
    append!(filled_points_p, p_circle[valid_points])
    
    # Add points along P_max
    if abs(p_max) < s_max
        push!(filled_points_q, q_at_pmax, -q_at_pmax)
        push!(filled_points_p, p_max, p_max)
    end
    
    plot!(plt, filled_points_q, filled_points_p, linecolor=:blue, linewidth=0, 
          fill=true, fillalpha=0.2, fillcolor=:lightblue)
    
    # Add origin point
    scatter!(plt, [0], [0], markersize=4, color=:black)
    
    return plt
end

"""
    plot_with_delta_limits(; p_min=-1.0, p_max=1.0, s_max=1.0, delta_max=60, eq=1.2, ut=1.0, xd=0.8, title="With Delta Limits")

Plot the cylindrical rotor model with delta limits (power angle constraints).
Q on x-axis, P on y-axis.
"""
function plot_with_delta_limits(; p_min=-1.0, p_max=1.0, s_max=1.0, delta_max=60, eq=1.2, ut=1.0, xd=0.8, title="With Delta Limits")
    # Convert delta_max from degrees to radians
    delta_max_rad = delta_max * π / 180
    delta_min_rad = -delta_max_rad
    
    # Create plot
    plt = plot(
        xlims=(-s_max*1.1, s_max*1.1),
        ylims=(p_min*1.1, p_max*1.1),
        xlabel="Q (pu)",
        ylabel="P (pu)",
        title=title,
        aspect_ratio=:equal,
        legend=false,
        grid=true
    )
    
    # Add y-axis line at x = 0
    plot!(plt, [0, 0], [p_min*1.1, p_max*1.1], linecolor=:black, linewidth=1, alpha=0.5)
    
    # Helper function to calculate P and Q for given delta
    function get_pq_from_delta(delta, eq_val)
        p = (eq_val * ut * sin(delta)) / xd
        q = (eq_val * ut * cos(delta) - ut^2) / xd
        return p, q
    end
    
    # Plot circular S_max constraint
    θ = range(0, 2π, length=100)
    q_circle = s_max .* cos.(θ)
    p_circle = s_max .* sin.(θ)
    
    # Apply P constraints to the circle
    valid_points = findall((p_circle .<= p_max) .& (p_circle .>= p_min))
    plot!(plt, q_circle[valid_points], p_circle[valid_points], linecolor=:blue, linewidth=2)
    
    # Calculate delta constraints using parametric approach to avoid zig-zag
    delta_range = range(-delta_max_rad, delta_max_rad, length=200)
    
    # Delta_max curve
    q_delta_max = []
    p_delta_max = []
    
    for delta in delta_range
        if abs(delta) <= delta_max_rad
            p, q = get_pq_from_delta(delta_max_rad, eq)
            
            # Check if point satisfies all constraints
            if p <= p_max && p >= p_min && p^2 + q^2 <= s_max^2
                push!(q_delta_max, q)
                push!(p_delta_max, p)
            end
        end
    end
    
    # Delta_min curve  
    q_delta_min = []
    p_delta_min = []
    
    for delta in delta_range
        if abs(delta) <= delta_max_rad
            p, q = get_pq_from_delta(delta_min_rad, eq)
            
            # Check if point satisfies all constraints
            if p <= p_max && p >= p_min && p^2 + q^2 <= s_max^2
                push!(q_delta_min, q)
                push!(p_delta_min, p)
            end
        end
    end
    
    # For the specific delta = +/- delta_max curves, we need to vary E_q
    eq_range = range(ut*0.8, eq*1.2, length=100)  # E_q should be higher than Ut
    
    # Delta = delta_max curve
    q_max_curve = []
    p_max_curve = []
    
    for eq_val in eq_range
        p, q = get_pq_from_delta(delta_max_rad, eq_val)
        
        # Check constraints
        if p <= p_max && p >= p_min && p^2 + q^2 <= s_max^2
            push!(q_max_curve, q)
            push!(p_max_curve, p)
        end
    end
    
    # Delta = delta_min curve
    q_min_curve = []
    p_min_curve = []
    
    for eq_val in eq_range
        p, q = get_pq_from_delta(delta_min_rad, eq_val)
        
        # Check constraints
        if p <= p_max && p >= p_min && p^2 + q^2 <= s_max^2
            push!(q_min_curve, q)
            push!(p_min_curve, p)
        end
    end
    
    # Plot delta constraints (now smooth curves)
    if !isempty(q_max_curve)
        plot!(plt, q_max_curve, p_max_curve, linecolor=:red, linewidth=2, label="δ = $delta_max °")
    end
    
    if !isempty(q_min_curve)
        plot!(plt, q_min_curve, p_min_curve, linecolor=:red, linewidth=2, label="δ = -$delta_max °")
    end
    
    # Add legend
    plot!(plt, legend=:topright)
    
    # Add origin point
    scatter!(plt, [0], [0], markersize=4, color=:black, label="")
    
    # Add annotation explaining the model
    annotate!(plt, s_max*0.5, p_min*0.9, text("Eq = $eq, Ut = $ut, Xd = $xd", :black, :left, 8))
    
    return plt
end

"""
    plot_full_constraints(; 
        p_min=-1.0, p_max=1.0, 
        q_min=-1.0, q_max=1.0,
        s_max=1.0, 
        delta_max=60, 
        eq_min=0.8, eq_max=1.2, 
        ut=1.0, xd=0.8, 
        title="Full Constraints Model")

Plot the complete generator capability model with all constraints:
- Apparent power limit (P^2 + Q^2 ≤ S_max^2)
- Active power limits (P_min ≤ P ≤ P_max)
- Reactive power limits (Q_min ≤ Q ≤ Q_max)
- Power angle limits (δ_min ≤ δ ≤ δ_max)
- Internal voltage limits (E_q,min ≤ E_q ≤ E_q,max)

Q on x-axis, P on y-axis.
"""
function plot_full_constraints(; 
    p_min=-1.0, p_max=1.0, 
    q_min=-1.0, q_max=1.0,
    s_max=1.0, 
    delta_max=60, 
    eq_min=0.8, eq_max=1.2, 
    ut=1.0, xd=0.8, 
    use_qmin_function=true,
    title="Full Constraints Model")
    
    # Ensure E_q values are higher than Ut for realistic generator operation
    eq_min = max(eq_min, ut * 1.05)  # E_q should be at least 5% higher than Ut
    eq_max = max(eq_max, ut * 1.1)   # E_q_max should be at least 10% higher than Ut
    
    # Convert delta_max from degrees to radians
    delta_max_rad = delta_max * π / 180
    delta_min_rad = -delta_max_rad
    
    # Create plot
    plt = plot(
        xlims=(q_min*1.2, q_max*1.2),
        ylims=(p_min*1.2, p_max*1.2),
        xlabel="Q (pu)",
        ylabel="P (pu)",
        title=title,
        aspect_ratio=:equal,
        grid=true
    )
    
    # Add y-axis line at x = 0
    plot!(plt, [0, 0], [p_min*1.2, p_max*1.2], linecolor=:black, linewidth=1, alpha=0.5)
    
    # Helper function to calculate P and Q for given delta and E_q
    function get_pq_from_delta_eq(delta, eq_val)
        p = (eq_val * ut * sin(delta)) / xd
        q = (eq_val * ut * cos(delta) - ut^2) / xd
        return p, q
    end
    
    # Plot circular S_max constraint
    θ = range(0, 2π, length=100)
    q_circle = s_max .* cos.(θ)
    p_circle = s_max .* sin.(θ)
    
    # Apply P and Q constraints to the circle
    valid_points = findall(
        (p_circle .<= p_max) .& 
        (p_circle .>= p_min) .& 
        (q_circle .<= q_max) .& 
        (q_circle .>= q_min)
    )
    
    plot!(plt, q_circle[valid_points], p_circle[valid_points], linecolor=:blue, linewidth=2, label="S ≤ $s_max")
    
    # Plot Q_min and Q_max vertical lines (only if they make sense)
    # For Q_min to be meaningful, we need to consider the physical generator limits
    # Typically Q_min is limited by field current or stability, not just a box constraint
    
    # Calculate realistic Q limits based on the generator equations
    # Q = (Eq*Ut*cos(delta) - Ut²)/Xd
    # For delta in [-delta_max, delta_max] and Eq in [eq_min, eq_max]
    
    q_theoretical_min = (eq_min * ut * cos(delta_max_rad) - ut^2) / xd
    q_theoretical_max = (eq_max * ut * cos(0) - ut^2) / xd  # cos(0) = 1 gives max Q
    
    # Use the more restrictive of theoretical and specified limits
    q_actual_min = max(q_min, q_theoretical_min)
    q_actual_max = min(q_max, q_theoretical_max)
    
    # Only plot Q limits if they are within the S constraint
    if q_actual_min > -s_max && q_actual_min < s_max
        p_at_qmin_upper = min(p_max, sqrt(s_max^2 - q_actual_min^2))
        p_at_qmin_lower = max(p_min, -sqrt(s_max^2 - q_actual_min^2))
        plot!(plt, [q_actual_min, q_actual_min], [p_at_qmin_lower, p_at_qmin_upper], 
              linecolor=:green, linewidth=2, label="Q_min")
    end
    
    if q_actual_max > -s_max && q_actual_max < s_max
        p_at_qmax_upper = min(p_max, sqrt(s_max^2 - q_actual_max^2))
        p_at_qmax_lower = max(p_min, -sqrt(s_max^2 - q_actual_max^2))
        plot!(plt, [q_actual_max, q_actual_max], [p_at_qmax_lower, p_at_qmax_upper], 
              linecolor=:green, linewidth=2, label="Q_max")
    end
    
    # Plot P_min and P_max horizontal lines
    if p_min > -s_max
        q_at_pmin_upper = min(q_actual_max, sqrt(s_max^2 - p_min^2))
        q_at_pmin_lower = max(q_actual_min, -sqrt(s_max^2 - p_min^2))
        plot!(plt, [q_at_pmin_lower, q_at_pmin_upper], [p_min, p_min], 
              linecolor=:purple, linewidth=2, label="P = $p_min")
    end
    
    if p_max < s_max
        q_at_pmax_upper = min(q_actual_max, sqrt(s_max^2 - p_max^2))
        q_at_pmax_lower = max(q_actual_min, -sqrt(s_max^2 - p_max^2))
        plot!(plt, [q_at_pmax_lower, q_at_pmax_upper], [p_max, p_max], 
              linecolor=:purple, linewidth=2, label="P = $p_max")
    end
    
    # Calculate and plot E_q limits curves using parametric approach
    delta_range = range(-delta_max_rad, delta_max_rad, length=200)
    
    # E_q,min curve
    q_eq_min = []
    p_eq_min = []
    
    for delta in delta_range
        p, q = get_pq_from_delta_eq(delta, eq_min)
        
        # Check if this point is within our other constraints
        if p <= p_max && p >= p_min && 
           q <= q_actual_max && q >= q_actual_min && 
           p^2 + q^2 <= s_max^2
            push!(q_eq_min, q)
            push!(p_eq_min, p)
        end
    end
    
    # E_q,max curve
    q_eq_max = []
    p_eq_max = []
    
    for delta in delta_range
        p, q = get_pq_from_delta_eq(delta, eq_max)
        
        # Check if this point is within our other constraints
        if p <= p_max && p >= p_min && 
           q <= q_actual_max && q >= q_actual_min && 
           p^2 + q^2 <= s_max^2
            push!(q_eq_max, q)
            push!(p_eq_max, p)
        end
    end
    
    # Plot E_q limits
    if !isempty(q_eq_min)
        plot!(plt, q_eq_min, p_eq_min, linecolor=:brown, linewidth=2, label="Eq = $eq_min")
    end
    
    if !isempty(q_eq_max)
        plot!(plt, q_eq_max, p_eq_max, linecolor=:orange, linewidth=2, label="Eq = $eq_max")
    end
    
    # Plot delta limits using parametric approach
    eq_range = range(eq_min, eq_max, length=100)
    
    # Delta = delta_max curve
    q_delta_max = []
    p_delta_max = []
    
    for eq_val in eq_range
        p, q = get_pq_from_delta_eq(delta_max_rad, eq_val)
        
        # Check constraints
        if p <= p_max && p >= p_min && 
           q <= q_actual_max && q >= q_actual_min && 
           p^2 + q^2 <= s_max^2
            push!(q_delta_max, q)
            push!(p_delta_max, p)
        end
    end
    
    # Delta = delta_min curve
    q_delta_min = []
    p_delta_min = []
    
    for eq_val in eq_range
        p, q = get_pq_from_delta_eq(delta_min_rad, eq_val)
        
        # Check constraints
        if p <= p_max && p >= p_min && 
           q <= q_actual_max && q >= q_actual_min && 
           p^2 + q^2 <= s_max^2
            push!(q_delta_min, q)
            push!(p_delta_min, p)
        end
    end
    
    # Plot delta constraints
    if !isempty(p_delta_max)
        plot!(plt, q_delta_max, p_delta_max, linecolor=:red, linewidth=2, label="δ = $delta_max°")
    end
    
    if !isempty(p_delta_min)
        plot!(plt, q_delta_min, p_delta_min, linecolor=:red, linewidth=2, label="δ = -$delta_max°")
    end
    
    # Add legend
    plot!(plt, legend=:topright)
    
    # Add origin point
    scatter!(plt, [0], [0], markersize=4, color=:black, label="")
    
    # Add annotation explaining the model parameters
    annotate!(plt, q_actual_min*0.8, p_min*0.9, 
              text("Eq: [$eq_min, $eq_max], Ut = $ut, Xd = $xd", :black, :left, 8))
    
    # Handle Q_min as function of P if requested
    if use_qmin_function
        # Calculate Q_min curve as function of P
        p_values = range(p_min, p_max, length=200)
        q_min_curve = []
        p_qmin_valid = []
        
        for p in p_values
            q_min_val = q_min_function(p, eq_min, ut, xd, delta_max_rad)
            if !isnan(q_min_val) && q_min_val >= q_min
                push!(q_min_curve, q_min_val)
                push!(p_qmin_valid, p)
            end
        end
        
        # Plot Q_min function
        if !isempty(q_min_curve)
            plot!(plt, q_min_curve, p_qmin_valid, linecolor=:red, linewidth=3, 
                  label="Q_min(P)")
        end
    else
        # Use constant Q_min (original approach)
        if q_min > -s_max && q_min < s_max
            p_at_qmin_upper = min(p_max, sqrt(s_max^2 - q_min^2))
            p_at_qmin_lower = max(p_min, -sqrt(s_max^2 - q_min^2))
            plot!(plt, [q_min, q_min], [p_at_qmin_lower, p_at_qmin_upper], 
                  linecolor=:green, linewidth=2, label="Q_min")
        end
    end
    
    return plt
end

"""
    compare_all_models(; kwargs...)

Compare all generator models in a single plot.
"""
function compare_all_models(; kwargs...)
    # Create individual plots
    simple_box = plot_simple_box_model(title="Simple Box")
    cylindrical_rotor = plot_cylindrical_rotor_model(title="Cylindrical Rotor")
    with_delta = plot_with_delta_limits(title="With Delta Limits")
    full = plot_full_constraints(title="Full Constraints")
    
    # Combine into a single figure
    combined = plot(simple_box, cylindrical_rotor, with_delta, full, layout=(2,2), size=(800, 800))
    plot!(combined, title="Generator Capability Curves Comparison")
    
    return combined
end

"""
    demonstrate_models()

Demonstrate the different generator capability models with examples.
"""
function demonstrate_models()
    # Example 1: Show individual models with default parameters
    simple_box = plot_simple_box_model()
    cylindrical_rotor = plot_cylindrical_rotor_model()
    with_delta = plot_with_delta_limits()
    full = plot_full_constraints()
    
    # Example 2: Compare all models
    comparison = compare_all_models()
    
    # Example 3: Full constraints model with custom parameters
    custom_model = plot_full_constraints(
        p_min=0.2, 
        p_max=1.0,
        q_min=-0.6,
        q_max=0.8,
        s_max=1.05,
        delta_max=45,
        eq_min=0.9,
        eq_max=1.3,
        ut=1.0,
        xd=0.7,
        title="Custom Full Constraints Model"
    )
    
    # Display plots
    display(simple_box)
    display(cylindrical_rotor)
    display(with_delta)
    display(full)
    display(comparison)
    display(custom_model)
    
    return nothing
end

"""
    q_min_function(p, eq_min, ut, xd, delta_max_rad)

Calculate Q_min as a function of P based on generator physical constraints.
This represents the minimum reactive power that can be produced at a given active power.
"""
function q_min_function(p, eq_min, ut, xd, delta_max_rad)
    # From P = (Eq*Ut/Xd)*sin(delta), we can find delta
    # Then Q = (Eq*Ut*cos(delta) - Ut²)/Xd
    
    # Calculate the minimum E_q needed for this P at maximum delta
    eq_needed_for_p = (p * xd) / (ut * sin(delta_max_rad))
    
    if eq_needed_for_p > eq_min
        # Use the required E_q at maximum delta
        delta = delta_max_rad
        eq_used = eq_needed_for_p
    else
        # Use minimum E_q and calculate corresponding delta
        eq_used = eq_min
        if abs(p) <= abs((eq_used * ut) / xd)
            delta = asin((p * xd) / (eq_used * ut))
        else
            return NaN  # P not achievable with this E_q
        end
    end
    
    # Calculate Q_min
    q_min = (eq_used * ut * cos(delta) - ut^2) / xd
    return q_min
end

"""
    plot_qmin_function(; p_min=-1.0, p_max=1.0, eq_min=1.05, ut=1.0, xd=0.8, delta_max=60, title="Q_min Function")

Plot Q_min as a function of P to show the curved constraint boundary.
"""
function plot_qmin_function(; p_min=-1.0, p_max=1.0, eq_min=1.05, ut=1.0, xd=0.8, delta_max=60, title="Q_min Function")
    delta_max_rad = delta_max * π / 180
    
    plt = plot(
        xlabel="Q (pu)",
        ylabel="P (pu)",
        title=title,
        aspect_ratio=:equal,
        grid=true,
        legend=:topright
    )
    
    # Add y-axis line at x = 0
    plot!(plt, [0, 0], [p_min*1.1, p_max*1.1], linecolor=:black, linewidth=1, alpha=0.5)
    
    # Calculate Q_min curve
    p_values = range(p_min, p_max, length=200)
    q_min_values = []
    p_valid = []
    
    for p in p_values
        q_min = q_min_function(p, eq_min, ut, xd, delta_max_rad)
        if !isnan(q_min)
            push!(q_min_values, q_min)
            push!(p_valid, p)
        end
    end
    
    # Plot the Q_min curve
    if !isempty(q_min_values)
        plot!(plt, q_min_values, p_valid, linecolor=:red, linewidth=3, label="Q_min(P)")
    end
    
    # Add some example operating points
    p_examples = [0.3, 0.6, 0.9]
    for p_ex in p_examples
        q_min_ex = q_min_function(p_ex, eq_min, ut, xd, delta_max_rad)
        if !isnan(q_min_ex)
            scatter!(plt, [q_min_ex], [p_ex], markersize=6, color=:blue, 
                    label=p_ex == p_examples[1] ? "Operating points" : "")
            annotate!(plt, q_min_ex + 0.1, p_ex, text("P=$p_ex", :blue, :left, 8))
        end
    end
    
    return plt
end

"""
    demonstrate_pv_to_pq_transition(; kwargs...)

Demonstrate the PV→PQ transition phenomenon when reactive power limits are hit.
"""
function demonstrate_pv_to_pq_transition(; 
    p_setpoint=0.8, 
    v_setpoint=1.0,
    eq_min=1.05, eq_max=1.4,
    ut=1.0, xd=0.8, 
    delta_max=60,
    title="PV → PQ Transition")
    
    delta_max_rad = delta_max * π / 180
    
    plt = plot(
        xlabel="Q (pu)",
        ylabel="P (pu)",
        title=title,
        aspect_ratio=:equal,
        grid=true,
        legend=:topright
    )
    
    # Add y-axis line at x = 0
    plot!(plt, [0, 0], [-0.2, 1.2], linecolor=:black, linewidth=1, alpha=0.5)
    
    # Calculate the capability curve
    delta_range = range(-delta_max_rad, delta_max_rad, length=200)
    
    # Upper boundary (Eq_max)
    q_upper = []
    p_upper = []
    
    for delta in delta_range
        p = (eq_max * ut * sin(delta)) / xd
        q = (eq_max * ut * cos(delta) - ut^2) / xd
        if p >= 0 && p <= 1.2  # Reasonable range
            push!(q_upper, q)
            push!(p_upper, p)
        end
    end
    
    # Lower boundary (Eq_min or Q_min function)
    q_lower = []
    p_lower = []
    
    for delta in delta_range
        p = (eq_min * ut * sin(delta)) / xd
        q = (eq_min * ut * cos(delta) - ut^2) / xd
        if p >= 0 && p <= 1.2
            push!(q_lower, q)
            push!(p_lower, p)
        end
    end
    
    # Plot capability boundaries
    if !isempty(q_upper)
        plot!(plt, q_upper, p_upper, linecolor=:blue, linewidth=2, label="Eq_max limit")
    end
    if !isempty(q_lower)
        plot!(plt, q_lower, p_lower, linecolor=:red, linewidth=2, label="Eq_min limit")
    end
    
    # Calculate Q limits for the given P setpoint
    # For PV mode: P is fixed, Q varies with system conditions
    p_fixed = p_setpoint
    
    # Calculate Q range for this P
    # From P = (Eq*Ut/Xd)*sin(delta), find delta for different Eq values
    eq_range = range(eq_min, eq_max, length=100)
    q_possible = []
    
    for eq in eq_range
        max_p_for_eq = (eq * ut) / xd
        if p_fixed <= max_p_for_eq
            delta = asin((p_fixed * xd) / (eq * ut))
            if abs(delta) <= delta_max_rad
                q = (eq * ut * cos(delta) - ut^2) / xd
                push!(q_possible, q)
            end
        end
    end
    
    if !isempty(q_possible)
        q_min_pv = minimum(q_possible)
        q_max_pv = maximum(q_possible)
        
        # Plot PV operating line
        plot!(plt, [q_min_pv, q_max_pv], [p_fixed, p_fixed], 
              linecolor=:green, linewidth=4, label="PV mode range")
        
        # Show transition points
        scatter!(plt, [q_min_pv], [p_fixed], markersize=8, color=:orange, 
                label="PV→PQ transition")
        scatter!(plt, [q_max_pv], [p_fixed], markersize=8, color=:orange, label="")
        
        # Annotations
        annotate!(plt, q_min_pv - 0.1, p_fixed + 0.05, 
                 text("Q_min hit\n(PV→PQ)", :orange, :right, 8))
        annotate!(plt, q_max_pv + 0.1, p_fixed + 0.05, 
                 text("Q_max hit\n(PV→PQ)", :orange, :left, 8))
        
        # Show a typical PV operating point
        q_normal = (q_min_pv + q_max_pv) / 2
        scatter!(plt, [q_normal], [p_fixed], markersize=6, color=:green, 
                label="Normal PV operation")
    end
    
    # Add explanatory text
    annotate!(plt, 0.5, 0.1, 
             text("PV mode: P fixed, Q varies\nPQ mode: Both P and Q fixed", :black, :center, 10))
    
    return plt
end

"""
    test_pv_pq_scenarios()

Test different scenarios showing PV→PQ transitions.
"""
function test_pv_pq_scenarios()
    # Scenario 1: Light load - generator can maintain voltage
    scenario1 = demonstrate_pv_to_pq_transition(
        p_setpoint=0.3,
        title="Scenario 1: Light Load (PV mode possible)"
    )
    
    # Scenario 2: Heavy load - Q_min limit hit
    scenario2 = demonstrate_pv_to_pq_transition(
        p_setpoint=0.8,
        title="Scenario 2: Heavy Load (Q_min limit hit)"
    )
    
    # Scenario 3: Very heavy load - Q_max limit hit
    scenario3 = demonstrate_pv_to_pq_transition(
        p_setpoint=0.5,
        eq_max=1.2,  # Lower max excitation
        title="Scenario 3: High Q demand (Q_max limit hit)"
    )
    
    # Combined plot
    combined = plot(scenario1, scenario2, scenario3, layout=(1,3), size=(1200, 400))
    
    return scenario1, scenario2, scenario3, combined
end

# Execute demonstration when the file is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    demonstrate_models()
end
