using LinearAlgebra, SparseArrays, Printf

x0 = zeros(2)
x0[1] = 1.5
x0[2] = 0.5

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