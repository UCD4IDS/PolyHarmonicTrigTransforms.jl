# LLST Laplace Local Sine Transform

# Applies LLST using given level list.

# Syntax:
#   coef = LLST(data, ll)

# Inputs:
#   data   - image to apply LLST to
#   ll   - level list

# Outputs:
#   coef - computed LLST coefficients

#   See also: ILLST.

module LLST

    include("dst.jl")
    include("solvelaplace.jl")

    using .DST
    using .SOLVELAPLACE
    
    export llst2d, illst2d

    # LLFT2D 2D Laplace Local Fourier Transform
    # 
    # Applies the LLFT transform to image (no blocking). Note that transform
    # is only applied to the interior of image. Boundary remains unchanged.
    # 
    # Syntax:
    #   transf = LLFT2D(data) 
    # 
    # Inputs:
    #   data     - image to be transformed
    # 
    # Outputs:
    #   transf - transformed image
    # 
    # See also: LLFT, ILLFT, SOLVELAPLACE.
    function llst2d(data::AbstractVecOrMat)
        (m, n) = size(data)

        trans = zeros(m, n)
        trans[:, 1] = data[:, 1]
        trans[:, end] = data[:, end]
        trans[end, :] = data[end, :]
        trans[1, :] = data[1, :]

        inIX1 = 2:m-1
        inIX2 = 2:n-1
        rnorm = sqrt(4 / (m - 1) / (n - 1))
        interior = data[inIX1, inIX2] .- solvelaplace(data)
        trans[inIX1, inIX2] = rnorm .* dst(dst(interior)')'
        return trans
    end

    # ILLFT2D 2D Inverse Laplace Local Fourier Transform
    # 
    # Applies the inverse LLFT transform to image (no blocking). Note that 
    # transform is only applied to the interior of image. Boundary remains 
    # unchanged.
    # 
    # Syntax:
    #   transf = ILLFT2D(im) 
    # 
    # Inputs:
    #   data     - image to be transformed
    # 
    # Outputs:
    #   transf - transformed image
    # 
    # See also: LLFT, ILLFT, SOLVELAPLACE.
    function illst2d(data::AbstractVecOrMat)
        (m, n) = size(data)

        trans = zeros(m, n)
        trans[:, 1] = data[:, 1]
        trans[:, end] = data[:, end]
        trans[end, :] = data[end, :]
        trans[1, :] = data[1, :]

        inIX1 = [2:1:m-1;]'
        inIX2 = [2:1:n-1;]'
        if ((m - 1) < 2)
            inIX1 = [2;]'
        end
        if ((n - 1) < 2)
            inIX2 = [2;]'
        end

        rnorm = sqrt(4 / (m - 1) / (n - 1))

        data = float(data)
        data[inIX1[:], inIX2[:]] = rnorm .* dst(dst(data[inIX1[:], inIX2[:]])')'
        trans[inIX1[:], inIX2[:]] = data[inIX1[:], inIX2[:]] .+ solvelaplace(trans)

        return trans
    end


end
