# SOLVELAPLACE solve the laplace equation.

# solves the laplace equation based on the data contained in 'im'. The size 
# of IM2 is (n-2)x(m-2) where n,m are the dimensions of 'im' (it returns 
# the "inside" of the solution as the borders are exactly the same as those 
# of 'im').

# Syntax:
#   im = SOLVELAPLACE(im)

# Inputs:
#   im - nxm matrix to solve laplace equation for

# Outputs:
#   im - (n-2)x(m-2) matrix laplacian solution. Boundary is of 'im' is
#        unchanged.

# See also SOLVELAPLACEOLD, LAPLACEEQ.

module SOLVELAPLACE

    include("dst.jl")
    include("idst.jl")
    using .DST
    using .IDST

    export solvelaplace

    function meshgrid(x::AbstractVecOrMat, y::AbstractVecOrMat)
        X = [x for _ in y, x in x]
        Y = [y for y in y, _ in x]
        X, Y
    end


    function solvelaplace(image::AbstractVecOrMat)
        if (isempty(image) == true)
            error("IM has to be a 2D array")
        end
        (m, n) = size(image)

        if (m == 2 & n == 2)
            return 0
        end

        d_x = 1 / (n - 1)
        d_y = 1 / (m - 1)
        #@info "d_x d_y" d_x d_y

        # [x,y] = meshgrid((1:n-2)*d_x,(1:m-2)*d_y);
        #x1 = (1:n-2)*d_x
        #y1 = (1:m-2)*d_y

        x, y = meshgrid((1:n-2) * d_x, (1:m-2) * d_y)
        #x = [x1 for _ in y1, x1 in x1]
        #y = [y1 for y1 in y1, _ in x1]
        #x = getindex.(Iterators.product(y1, x1), 2)
        #y = getindex.(Iterators.product(y1, x1), 1)

        #@info "x y" x y

        b1 = image[2:end-1, 1]
        b2 = image[end, 2:end-1]'
        b3 = image[2:end-1, end]
        b4 = image[1, 2:end-1]'

        # Remove a function s.t. the laplacian of the borders

        # Remove the 4 corners,
        # z = axy+bx+cy+d with x,y in [0,1]x[0,1]

        d = image[1, 1]
        b = image[1, n] - d
        c = image[m, 1] - d
        a = image[m, n] - b - c - d

        #lapl = a .* x .* y .+ b .* x .+ c .* y .+ d;
        lapl = a * x .* y .+ b .* x .+ c .* y .+ d

        b1 = b1 - (c * [1:m-2;] * d_y) .- d
        b4 = b4 - (b * [1:n-2;] * d_x)' .- d
        b3 = b3 .- b .- ((a + c) * [1:m-2;] * d_y) .- d
        b2 = b2 .- c .- ((a + b) * [1:n-2;] * d_x)' .- d

        # Now for the actual "solution" of the poisson equation
        db1 = dst(b1)
        db2 = dst(b2, 2)
        db3 = dst(b3)
        db4 = dst(b4, 2)

        #@info "dst" db1 db2 db3 db4

        lambdan = pi * [1:n-2;]'
        lambdam = pi * [1:m-2;]

        # if x>=20, then (sinh(x)-exp(x)/2)/exp(x) is below the double precision.
        # we will use this threshold to switch from the sinh value to the
        # exponential approximation

        #@info "lambda" lambdan lambdam

        u = repeat(lambdam, 1, n - 2)
        v = repeat(lambdan, m - 2, 1)

        #@info "h" x u v

        # h1 = sinh(repmat(lambdam,1,n-2).*(1-x))./sinh(repmat(lambdam,1,n-2));
        #h1 = exp(-u.*x) .* ( 1 .- exp.(-2 * u .*(1 .- x)))./(1 .- exp.(-2*u));
        h1 = exp.(-u .* x) .* (1 .- exp.(-2 .* u .* (1 .- x))) ./ (1 .- exp.(-2 .* u))
        # h3 = sinh(repmat(lambdam,1,n-2).*x)./sinh(repmat(lambdam,1,n-2));
        h3 = exp.(u .* (x .- 1)) .* (1 .- exp.(-2 * u .* x)) ./ (1 .- exp.(-2 * u))
        # h2 = sinh(repmat(lambdan,m-2,1).*y)./sinh(repmat(lambdan,m-2,1));
        h2 = exp.(v .* (y .- 1)) .* (1 .- exp.(-2 * v .* y)) ./ (1 .- exp.(-2 * v))
        # h4 = sinh(repmat(lambdan,m-2,1).*(1-y))./sinh(repmat(lambdan,m-2,1));
        h4 = exp.(-v .* y) .* (1 .- exp.(-2 * v .* (1 .- y))) ./ (1 .- exp.(-2 * v))

        #@info "repeat" db1 db2 db3 db4 
        #@info "repeat" m n h1 h2 h3 h4
        h1m = repeat(db1, 1, n - 2) .* h1
        h3m = repeat(db3, 1, n - 2) .* h3
        h2m = repeat(db2, m - 2, 1) .* h2
        h4m = repeat(db4, m - 2, 1) .* h4

        #@info "last" lapl h1m h2m h3m h4m

        hi1 = h1m + h3m
        i1 = idst(hi1)
        i2 = idst(h2m' + h4m')'
        lapl = lapl + i1 + i2
        #lapl = lapl + idst(h1m+h3m) + idst(h2m'+h4m')';
        return lapl
    end
end;