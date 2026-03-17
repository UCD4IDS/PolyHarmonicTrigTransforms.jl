# IDST Inverse Discrete Sine Transform

# Computes the inverse type I discrete sine transform of 'y'.  If 'n' is 
# given, then 'y' is padded or trimmed to length 'n' before computing the 
# transform. If 'y' is a matrix, compute the transform along the columns of 
# the the matrix.

# Syntax:
#   x = IDST(y)
#   x = IDST(y, n)

# Inputs:
#   x - vector to be transformed
#   n - length of tranform to perform

# Outputs:
#   y - tranformed vector

# Author:
#   Paul Kienzle <pkienzle@users.sf.net> (2006)
#   This program is granted to the public domain.

# See also: DST
# 
module IDST
    include("dst.jl")
    using .DST

    # Lazily load FFTW when needed by IDST helpers.
    function _ensure_fftw()
        if !isdefined(@__MODULE__, :FFTW)
            Core.eval(Main, :(using FFTW))
            Core.eval(@__MODULE__, :(const FFTW = Main.FFTW))
        end
        return nothing
    end

    export idst, idst!, plan_idst, idst_old
    """
    idst(y; dims=1)

    Inverse Type-I Discrete Sine Transform (IDST) of `y` along `dims`.

    Arguments
    - `y`: an `AbstractArray` of DST coefficients.
    - `dims`: dimension to transform (default `1`).

    Returns
    - reconstructed array in the original domain.
    """
    function idst(y::AbstractArray, dims=1)
        
        N = size(y, dims)
        x = dst(y / (2*(N+1)), dims)
        #x = dst(y * 2/ (N+1), dims)
        #@info "idst" x
        return x * 4
    end

    """
    plan_idst(y; dims=1, flags=FFTW.MEASURE)

    Create an FFTW r2r plan suitable for the IDST (uses the DST plan internally).
    """
    function plan_idst(y::AbstractArray; dims=1, flags=FFTW.MEASURE)
        _ensure_fftw()
        return FFTW.r2r_plan(y, FFTW.RODFT00, dims; flags=flags)
    end

    """
    idst!(y; dims=1)

    In-place inverse DST: overwrite `y` with the reconstructed data when
    possible. Falls back to computing out-of-place and copying the result.
    Uses an FFTW r2r plan for the in-place computation when available. The result is scaled by `1/(N+1)` to account for the normalization of the DST. Note that the input `y` should contain the DST coefficients, and after this function, it will contain the reconstructed data in the original domain.
    """
    function idst!(y::AbstractArray; dims=1)
        _ensure_fftw()
        try
            # In-place inverse DST: perform r2r in-place then scale by 1/(N+1)
            FFTW.r2r!(y, FFTW.RODFT00, dims)
            y ./= (size(y, dims) + 1)
            return y
        catch _
            tmp = idst(y, dims)
            y .= tmp
            return y
        end
    end

    function idst_old(y::AbstractArray, n=nothing, dims=1)
        if n === nothing
            n = size(y, 1)
            if n == 1
                n = size(y, 2)
            end
        end

        x = dst_old(y, n)
        x = x .* 2 / (n + 1)
        #@info "idst" x
        return x
    end
    # !test
    # ! x = log(gausswin(32));
    # ! assert(x, idst(dst(x)), 100*eps)
end;

