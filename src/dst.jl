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
    using AbstractFFTs, FFTW
    
    export dst

    function dst(a::AbstractArray, n=nothing)
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
        yy = fft(y, (1,))
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