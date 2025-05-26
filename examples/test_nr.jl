"""
Newton-Raphson Test for System of Equations
Solving:
x₁² - x₂ - 1 = 0
-x₁ + x₂² - 1 = 0
"""

using LinearAlgebra, SparseArrays, Printf

function test_system_newton_raphson()
    println("Testing Newton-Raphson for system of equations:")
    println("x₁² - x₂ - 1 = 0")
    println("-x₁ + x₂² - 1 = 0")
    
    x0 = [1.5, 0.5]  # Initial guess
    
    function create_jacobian(x)
        J = zeros(2, 2)
        J[1, 1] = 2 * x[1]
        J[1, 2] = -1
        J[2, 1] = -1
        J[2, 2] = 2 * x[2]
        return J
    end

    function create_residual(x)
        r = zeros(2)
        r[1] = x[1]^2 - x[2] - 1
        r[2] = -x[1] + x[2]^2 - 1
        return r
    end

    function newton_raphson(x0, tol=1e-6, max_iter=100)
        x = copy(x0)
        println("Starting Newton-Raphson iteration...")
        println("Initial guess: x = $x")
        
        for i in 1:max_iter
            println("Iteration $i:")
            J = create_jacobian(x)
            println("  Jacobian:")
            display(J)
            r = create_residual(x)
            println("  Residual: $r")
            
            if norm(r) < tol
                println("✓ Converged in $i iterations!")
                println("Solution: x₁ = $(x[1]), x₂ = $(x[2])")
                return x, true, i
            end
            
            dx = -J \ r
            x += dx
            println("  New x: $x")
            println()
        end
        
        println("✗ Failed to converge in $max_iter iterations")
        return x, false, max_iter
    end

    return newton_raphson(x0)
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    test_system_newton_raphson()
end
