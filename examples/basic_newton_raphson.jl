"""
Basic Newton-Raphson Examples
Demonstrates Newton-Raphson method for different equation types
"""

using LinearAlgebra, SparseArrays, Printf

# Example 1: Polynomial equation x³ - 2x - 5 = 0
function example_polynomial()
    println("Example 1: Solving x³ - 2x - 5 = 0")
    
    x0 = [-1000.0]
    
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
        for i in 1:max_iter
            J = create_jacobian(x)
            r = create_residual(x)
            if norm(r) < tol
                println("Converged in $i iterations: x = $(x[1])")
                return x, true, i
            end
            dx = -J \ r
            x += dx
        end
        println("Failed to converge")
        return x, false, max_iter
    end
    
    return newton_raphson(x0)
end

# Example 2: System of equations
function example_system()
    println("Example 2: System of equations")
    println("x₁² - x₂ - 1 = 0")
    println("-x₁ + x₂² - 1 = 0")
    
    x0 = [1.5, 0.5]
    
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
        for i in 1:max_iter
            J = create_jacobian(x)
            r = create_residual(x)
            if norm(r) < tol
                println("Converged in $i iterations: x₁ = $(x[1]), x₂ = $(x[2])")
                return x, true, i
            end
            dx = -J \ r
            x += dx
        end
        println("Failed to converge")
        return x, false, max_iter
    end
    
    return newton_raphson(x0)
end

# Example 3: Divergent case
function example_divergent()
    println("Example 3: Potentially divergent case - tan(x) - x = 0")
    
    using Plots
    
    x0 = [10.0]
    
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
        for i in 1:max_iter
            J = create_jacobian(x)
            r = create_residual(x)
            if norm(r) < tol
                println("Converged in $i iterations: x = $(x[1])")
                return x, true, i
            end
            dx = -J \ r
            x += dx
            println("Iteration $i: x = $(x[1])")
        end
        println("Failed to converge within $max_iter iterations")
        return x, false, max_iter
    end
    
    # Plot the function
    x_plot = range(-10, 10, length=1000)
    y_plot = tan.(x_plot) - x_plot
    p = plot(x_plot, y_plot, label="f(x) = tan(x) - x", 
             xlabel="x", ylabel="f(x)", title="f(x) = tan(x) - x", 
             lw=2, ylims=(-10, 10))
    savefig(p, "tan_function.png")
    
    return newton_raphson(x0)
end

# Run all examples
function main()
    println("="^50)
    println("BASIC NEWTON-RAPHSON EXAMPLES")
    println("="^50)
    
    example_polynomial()
    println()
    example_system()
    println()
    example_divergent()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
