"""
Newton-Raphson Test for Potentially Divergent Case
Solving f(x) = tan(x) - x = 0
"""

using LinearAlgebra, SparseArrays, Printf
using Plots

function test_divergent_newton_raphson()
    println("Testing Newton-Raphson for potentially divergent case: f(x) = tan(x) - x = 0")
    
    x0 = [10.0]  # Initial guess
    
    function create_jacobian(x)
        J = zeros(1, 1)
        J[1, 1] = sec(x[1])^2 - 1
        return J
    end

    function create_residual(x)
        r = zeros(1)
        r[1] = tan(x[1]) - x[1]
        return r
    end

    function newton_raphson(x0, tol=1e-6, max_iter=20)
        x = copy(x0)
        println("Starting Newton-Raphson iteration...")
        println("Initial guess: x = $(x[1])")
        println("Note: This function often diverges due to discontinuities in tan(x)")
        
        for i in 1:max_iter
            println("Iteration $i:")
            J = create_jacobian(x)
            println("  Jacobian: $J")
            r = create_residual(x)
            println("  Residual: $r")
            
            if norm(r) < tol
                println("✓ Converged in $i iterations!")
                println("Solution: x = $(x[1])")
                return x, true, i
            end
            
            if abs(J[1,1]) < 1e-12
                println("✗ Jacobian is nearly singular, stopping iteration")
                return x, false, i
            end
            
            dx = -J \ r
            x += dx
            println("  New x: $(x[1])")
            
            # Check for unreasonable values
            if abs(x[1]) > 1000
                println("✗ Solution diverged (|x| > 1000)")
                return x, false, i
            end
            println()
        end
        
        println("✗ Failed to converge in $max_iter iterations")
        return x, false, max_iter
    end

    # Plot the function
    println("\nPlotting f(x) = tan(x) - x...")
    x_plot = range(-10, 10, length=1000)
    y_plot = tan.(x_plot) - x_plot
    
    p = plot(x_plot, y_plot, 
             label="f(x) = tan(x) - x", 
             xlabel="x", 
             ylabel="f(x)", 
             title="f(x) = tan(x) - x", 
             lw=2, 
             ylims=(-10, 10),
             size=(800, 600))
    
    # Add horizontal line at y=0
    hline!([0], color=:red, linestyle=:dash, label="y = 0")
    
    # Save plot
    savefig(p, "tan_function.png")
    println("Function plot saved as 'tan_function.png'")
    
    return newton_raphson(x0)
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    result = test_divergent_newton_raphson()
    println("\nTest completed.")
end
