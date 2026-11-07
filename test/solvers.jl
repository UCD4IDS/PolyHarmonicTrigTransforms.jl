# Accuracy tests for the rectangle solvers `solvelaplace` (Dirichlet) and
# `solveneumann` (Neumann, constant right-hand side), against exact solutions.

using Test

@testset "solvelaplace (Dirichlet)" begin
    # A harmonic bilinear function is reproduced exactly (corner term only).
    m, n = 17, 33
    X = [(i - 1) / (n - 1) for j in 1:m, i in 1:n]; Y = [(j - 1) / (m - 1) for j in 1:m, i in 1:n]
    g = 1 .+ 2X .- 3Y .+ 5X .* Y
    @test maximum(abs.(solvelaplace(g) .- g)) < 1e-12

    # Square grid: e^x cos y on the unit square, second-order accurate.
    for (N, tol) in ((33, 2e-6), (129, 1e-7))
        X = [(i - 1) / (N - 1) for j in 1:N, i in 1:N]; Y = [(j - 1) / (N - 1) for j in 1:N, i in 1:N]
        g = exp.(X) .* cos.(Y)
        @test maximum(abs.(solvelaplace(g) .- g)) < tol
    end

    # Non-square grid with square pixels: the rectangle has aspect ratio
    # r = (n-1)/(m-1) = 2, and exp(r x) cos(y) is harmonic in physical coordinates.
    m, n = 65, 129; r = (n - 1) / (m - 1)
    X = [(i - 1) / (n - 1) for j in 1:m, i in 1:n]; Y = [(j - 1) / (m - 1) for j in 1:m, i in 1:n]
    g = exp.(r .* X) .* cos.(Y)
    @test maximum(abs.(solvelaplace(g) .- g)) < 1e-5
    # The stretched option reproduces the pre-0.1.7 behavior: harmonic in (X, Y).
    g2 = exp.(X) .* cos.(Y)
    @test maximum(abs.(solvelaplace(g2; stretched=true) .- g2)) < 1e-5
    # Square grids are unaffected by the option.
    s = randn(9, 9)
    @test solvelaplace(s) ≈ solvelaplace(s; stretched=true)
    # Degenerate sizes have no interior.
    @test solvelaplace([1.0 2.0; 3.0 4.0]) == [1.0 2.0; 3.0 4.0]
end

@testset "solveneumann (Neumann)" begin
    # Exact solution w = e^x cos(y) + 0.3 (x^2 + y^2), with Laplacian K = 1.2,
    # on (0,a) x (0,b) with a = 1 + sqrt(2), b = 1.
    a, b = 1 + sqrt(2), 1.0
    w(x, y)  = exp(x) * cos(y) + 0.3 * (x^2 + y^2)
    wx(x, y) = exp(x) * cos(y) + 0.6x
    wy(x, y) = -exp(x) * sin(y) + 0.6y
    # Measured (max over the grid / middle half of the grid):
    #   (32, 64): 2.39e-3 / 3.41e-4;  (128, 256): 2.01e-4 / 2.16e-5.
    # Second order in the interior; the maximum sits at the corner (a, b).
    errs = Float64[]; ierrs = Float64[]
    for (m, n) in ((32, 64), (128, 256))
        x = ((0:n-1) .+ 0.5) .* (a / n); y = ((0:m-1) .+ 0.5) .* (b / m)
        u, K = solveneumann(-wx.(0, y), wx.(a, y), -wy.(x, 0), wy.(x, b); a=a, b=b)
        @test size(u) == (m, n)
        @test isapprox(K, 1.2; atol=1e-3)
        W = [w(xi, yj) for yj in y, xi in x]
        d = u .- W
        d .-= sum(d) / length(d)                                # equal up to a constant
        push!(errs, maximum(abs.(d)))
        push!(ierrs, maximum(abs.(d[m÷4:3m÷4, n÷4:3n÷4])))
    end
    @test errs[2] < 3e-4
    @test ierrs[2] < 3e-5
    @test ierrs[1] / ierrs[2] > 14                            # second order: 4^2 = 16
    @test errs[1] / errs[2] > 10

    # Compatible data with zero means on every edge give K = 0 and a harmonic u.
    m, n = 16, 16
    y = ((0:m-1) .+ 0.5) ./ m; x = ((0:n-1) .+ 0.5) ./ n
    u, K = solveneumann(cos.(pi .* y), zeros(m), zeros(n), zeros(n))
    @test abs(K) < 1e-14
    # Exact solution: u = cos(pi y) cosh(pi (1 - x)) / (pi sinh(pi)).
    @test maximum(abs.(u .- [cos(pi * yj) * cosh(pi * (1 - xi)) / (pi * sinh(pi)) for yj in y, xi in x])) < 1e-12
end
