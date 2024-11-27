module ILLST
    # ILLST Inverse Laplace Local Sine Transform

    # Applies inverse LLST using given level list.

    # Syntax:
    #   image = ILLST(coef, ll)

    # Inputs:
    #   coef - matrix of coefficients computed with LLST
    #   ll   - level list

    # Outputs:
    #   image   - reconstructed image

    # See also: LLST.
    include("llst.jl")
    using .LLST
    
    export illst
    
    function illst(image::AbstractMatrix, ll::AbstractMatrix)
        image = llst(image, ll, true)
    end
end