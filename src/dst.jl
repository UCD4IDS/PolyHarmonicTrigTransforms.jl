# For Intel might need to run for MKL 
# FFTW.set_provider!("mkl")
# DST Discrete Sine Transform.
# 
# Computes the type I discrete sine transform of 'x'.  If 'n' is given, 
# then 'x' is padded or trimmed to length 'n' before computing the 
# transform. If 'x' is a matrix, compute the transform along the columns of 
# the matrix.
# 
# The discrete sine transform X of x can be defined as follows:
#        N
# X[k] = sum x[n] sin (pi n k / (N+1)),  k = 1, ..., N
#        n=1
# 
# Syntax:
#   y = DST(x)
#   y = DST(x, n)
# 
# Inputs:
#   x - vector to be transformed
#   n - length of tranform to perform
# 
# Outputs:  
#   y - tranformed vector
# 
# Author:
#   Jason Hee Yoon (2022)
#   Shuchen Ye (2022)

# See also: IDST
module DST
    function _ensure_fftw()
        if !isdefined(@__MODULE__, :FFTW)
            Core.eval(Main, :(using FFTW))
            Core.eval(@__MODULE__, :(const FFTW = Main.FFTW))
        end
        return nothing
    end
    
    export dst, dst!, dst_u, plan_dst, dst_old

    """
        dst(x; dims=1)

    Compute the Type-I Discrete Sine Transform (DST) of `x` along `dims`.

    Arguments
    - `x`: an `AbstractArray` (vector or matrix).
    - `dims`: dimension to transform (default `1`).

    Returns
    - transformed array of the same shape as `x`.

    Example
    ```julia
    y = dst([1.0,2.0,3.0])
    ```
    """
    function dst(x::AbstractArray, dims=1)
        _ensure_fftw()
        N = size(x, dims)
        return FFTW.r2r(x, FFTW.RODFT00, dims) / 2
    end

    function dst_u(x::AbstractArray, dims=1)
        _ensure_fftw()
        N = size(x, dims)
        new_x = FFTW.r2r(x, FFTW.RODFT00, dims)

        return new_x / (sqrt(2*N)) #unitory value
    end

    """
        plan_dst(x; dims=1, flags=FFTW.MEASURE)

    Create an FFTW r2r plan for the DST (RODFT00) on `x`. The returned plan
    can be executed with `plan * x`.
    """
    function plan_dst(x::AbstractArray; dims=1, flags=FFTW.MEASURE)
        _ensure_fftw()
        return FFTW.r2r_plan(x, FFTW.RODFT00, dims; flags=flags)
    end

    """
        dst!(x; dims=1)

    In-place DST that overwrites `x` with its transform when possible. Falls
    back to a safe out-of-place compute and copy if the in-place routine is
    not available on the current platform.
    """
    function dst!(x::AbstractArray; dims=1)
        _ensure_fftw()
        try
            # Try to use FFTW in-place r2r if available
            FFTW.r2r!(x, FFTW.RODFT00, dims)
            x .*= 0.5
            return x
        catch _
            # Fallback: compute out-of-place and copy into x
            tmp = FFTW.r2r(x, FFTW.RODFT00, dims)
            tmp .*= 0.5
            x .= tmp
            return x
        end
    end

    function dst_old(a::AbstractArray, n=nothing)
        _ensure_fftw()
        if minimum(size(a)) == 1
            if size(a, 2) > 1
                do_trans = 1
            else
                do_trans = 0
            end
            a = a[:]
        else
            do_trans = 0
        end
        if n === nothing
            n = size(a, 1)
        end
        m = size(a, 2)

        # Pad or truncate a if necessary
        if size(a, 1) < n
            aa = zeros(n, m)
            aa[1:size(a, 1), :] = a
        else
            aa = a[1:n, :]
        end

        y = zeros(2 * (n + 1), m)
        y[2:n+1, :] = aa
        y[n+3:2*(n+1), :] = -reverse(aa, dims=1)
        yy = FFTW.fft(y, (1,))
        b = yy[2:n+1, :] / (-2 * sqrt(-1 + 0im))
        

        if isreal(a)
            b = real(b)
            #b[b.<0] .=0
        end
        if do_trans == 1
            b = b'
        end
        return b
    end
end;