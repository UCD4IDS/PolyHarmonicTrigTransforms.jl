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

    export idst
    function idst(y::AbstractArray, n=nothing)
        if n === nothing
            n = size(y, 1)
            if n == 1
                n = size(y, 2)
            end
        end

        x = dst(y, n)
        x = x .* 2 / (n + 1)
        #@info "idst" x
        return x
    end
    # !test
    # ! x = log(gausswin(32));
    # ! assert(x, idst(dst(x)), 100*eps)
end;

