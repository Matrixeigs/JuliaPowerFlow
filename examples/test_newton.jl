"""
Basic Newton-Raphson Test for Polynomial Equation
Solving x³ - 2x - 5 = 0
"""

using LinearAlgebra, SparseArrays, Printf

function test_polynomial_newton_raphson()
    println("Testing Newton-Raphson for polynomial equation: x³ - 2x - 5 = 0")
    
    x0 = [-1000.0]  # Initial guess
    
    function create_jacobian(x)
        J = zeros(1, 1)
        J[1, 1] = 3 * x[1]^2 - 2
        return J
    end

    function create_residual(x)
        r = zeros(1)
        r[1] = x[1]^3 - 2*x[1] - 5
        return r
    end

    function newton_raphson(x0, tol=1e-6, max_iter=100)
        x = copy(x0)
        println("Starting Newton-Raphson iteration...")
        println("Initial guess: x = $(x[1])")
        
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
            
            dx = -J \ r
            x += dx
            println("  New x: $(x[1])")
            println()
        end
        
        println("✗ Failed to converge in $max_iter iterations")
        return x, false, max_iter
    end

    return newton_raphson(x0)
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    test_polynomial_newton_raphson()
end
