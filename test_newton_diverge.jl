using LinearAlgebra, SparseArrays, Printf
# f(x) = tan(x) - x

x0 = zeros(1)
x0[1] = 10.0

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

function newton_raphson(x0, tol=1e-6, max_iter=1000)
    x = x0
    for i in 1:max_iter
        @show i
        J = create_jacobian(x)
        @show J
        r = create_residual(x)
        @show r
        dx = -J \ r
        x += dx
        @show x
        if norm(dx) < tol
            break
        end
    end
    return x
end

newton_raphson(x0)