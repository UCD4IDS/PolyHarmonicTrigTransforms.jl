# SOLVENEUMANN solve the Neumann problem for Poisson's equation with a constant
# right-hand side on a rectangle.
#
# Given the outward normal derivative h on the four edges of the rectangle
# (0,a) x (0,b), it computes the solution u of
#
#     Laplacian(u) = K  in the rectangle,     du/dnu = h  on its boundary,
#
# where K = (1/(a*b)) * (integral of h over the boundary) is the only constant for
# which this problem is solvable. This is the u-component of the Neumann
# decomposition f = u + v of Yamatani and Saito (PHLCT), with h = df/dnu: the residual
# v = f - u then satisfies the homogeneous Neumann condition.
#
# Construction: u = q + sum of four harmonic edge series, where the quadratic
#     q(x,y) = alpha*x^2 + beta*y^2 + gamma*x + delta*y,   with Laplacian(q) = K,
# absorbs the mean of h on each edge, and each edge series
#     sum_k c_k cos(pi*k*s/L) cosh(...)/(kappa*sinh(...))
# matches the remaining (zero-mean) data on its edge with zero normal derivative on
# the other three. The edge data are sampled at the midpoints of the grid cells,
# so the coefficients c_k are DCT-II coefficients, and u is returned on the
# midpoint grid of the rectangle, ready for the DCT-II of v = f - u.
#
# Author: Naoki Saito's group, 2026.

module SOLVENEUMANN

    using FFTW

    export solveneumann

    # cosh(kappa*s)/(kappa*sinh(kappa*L)) without overflow, for kappa > 0 and 0 <= s <= L.
    _coshratio(kappa, s, L) =
        (exp(kappa * (s - L)) + exp(-kappa * (s + L))) / (kappa * (1 - exp(-2 * kappa * L)))

    # Coefficients c_k, k = 1, ..., N-1, of the cosine series of zero-mean data
    # sampled at the midpoints (j + 1/2) L/N: c_k = (2/N) sum_j h_j cos(pi k (j+1/2)/N).
    _cosinecoefs(h) = FFTW.r2r(h, FFTW.REDFT10)[2:end] ./ length(h)

    # Sum_k c_k cos(pi k (j+1/2)/N) for j = 0, ..., N-1, i.e., a DCT-III with zero DC term,
    # applied along dimension `dims` of the coefficient array C (k = 1, ..., N-1).
    function _cosinesynth(C::AbstractMatrix, dims::Int)
        Z = dims == 1 ? vcat(zeros(1, size(C, 2)), C) : hcat(zeros(size(C, 1), 1), C)
        return FFTW.r2r(Z ./ 2, FFTW.REDFT01, dims)
    end

    """
        solveneumann(hW, hE, hS, hN; a=1.0, b=1.0)

    Solve `Δu = K` in the rectangle `(0,a) × (0,b)` with outward normal derivative
    `∂u/∂ν = h` on its boundary, where `K = (∫ h dσ)/(ab)` is the unique constant that
    makes the problem solvable. In the PHLCT decomposition `f = u + v`, pass
    `h = ∂f/∂ν`; then `v = f - u` satisfies the homogeneous Neumann condition.

    Arguments
    - `hW`, `hE`: outward normal derivatives on the edges `x = 0` and `x = a`
      (so `hW = -∂u/∂x`, `hE = ∂u/∂x`), sampled at `y_j = (j + 1/2) b/m`, `j = 0, …, m-1`.
    - `hS`, `hN`: outward normal derivatives on the edges `y = 0` and `y = b`
      (so `hS = -∂u/∂y`, `hN = ∂u/∂y`), sampled at `x_i = (i + 1/2) a/n`, `i = 0, …, n-1`.
    - `a`, `b`: width and height of the rectangle.

    Returns `(u, K)`, where `u` is the m×n array of the solution at the cell midpoints
    `(x_i, y_j)`, with rows indexed by `j` (the `y` direction) and columns by `i`.
    The solution is unique up to an additive constant; the one returned has the mean
    of its four edge series equal to zero.
    """
    function solveneumann(hW::AbstractVector, hE::AbstractVector,
                          hS::AbstractVector, hN::AbstractVector;
                          a::Real=1.0, b::Real=1.0)
        m, n = length(hW), length(hS)
        length(hE) == m || error("hW and hE must have the same length")
        length(hN) == n || error("hS and hN must have the same length")
        (m >= 1 && n >= 1) || error("edge data must be nonempty")
        x = ((0:n-1) .+ 0.5) .* (a / n)
        y = ((0:m-1) .+ 0.5) .* (b / m)

        # Edge means (midpoint rule) and the quadratic q that absorbs them.
        mW, mE, mS, mN = sum(hW) / m, sum(hE) / m, sum(hS) / n, sum(hN) / n
        alpha, gamma = (mE + mW) / (2a), -mW
        beta,  delta = (mN + mS) / (2b), -mS
        K = 2alpha + 2beta              # = (b*(mW + mE) + a*(mS + mN)) / (a*b)
        u = [alpha * xi^2 + beta * yj^2 + gamma * xi + delta * yj for yj in y, xi in x]

        # West and east edges: cosine series in y, decaying across x.
        if m > 1
            cW, cE = _cosinecoefs(hW .- mW), _cosinecoefs(hE .- mE)
            kap = pi .* (1:m-1) ./ b
            C = [cW[k] * _coshratio(kap[k], a - xi, a) + cE[k] * _coshratio(kap[k], xi, a)
                 for k in 1:m-1, xi in x]
            u .+= _cosinesynth(C, 1)
        end
        # South and north edges: cosine series in x, decaying across y.
        if n > 1
            cS, cN = _cosinecoefs(hS .- mS), _cosinecoefs(hN .- mN)
            kap = pi .* (1:n-1) ./ a
            C = [cS[k] * _coshratio(kap[k], b - yj, b) + cN[k] * _coshratio(kap[k], yj, b)
                 for yj in y, k in 1:n-1]
            u .+= _cosinesynth(C, 2)
        end
        return u, K
    end
end
