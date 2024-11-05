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
    
    export llst, llst2D, illst2Ds
        
    function llst(data::AbstractVecOrMat, ll::AbstractArray, inverse::Bool=false)
        if (inverse != false)
            inverse = true
        end

        if (inverse)  # Inverse transform

            # Now we find the boundaries and process them.
            data = gridBlkBdryTrans(data, ll, process_border, inverse)

            # Recover interior of regions
            data = gridBlkTrans(data, ll, illst2D)

        else  # Forward transform
            # First we transform the interior of each region leaving boundaries
            # untouched.
            data = gridBlkTrans(data, ll, llst2D)

            # Now we find the boundaries and process them.
            data = gridBlkBdryTrans(data, ll, process_border, inverse)
        end
    end

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
    function llst2D(data::AbstractVecOrMat)
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
    function illst2D(data::AbstractVecOrMat)
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

    # Function for forward transform of boundary
    function fbdrytrans(bdry::AbstractVecOrMat)
        bdry = remove_line(bdry)
        bdry = sqrt(2 / (length(bdry) + 1)) * dst(bdry)
    end

    # Function for inverse transform of boundary
    function ibdrytrans(bdry::AbstractVecOrMat)
        x = sqrt(2 / (length(bdry) - 1)) * dst(bdry[2:end-1]')
        bdry = [bdry[1] x bdry[end]]
        bdry = addb_line(bdry)
    end

    function process_border(border::AbstractVecOrMat, ix::AbstractVecOrMat, fn::Function)
        for i = 1:length(ix)-1
            b = ix[i]
            e = ix[i+1]
            if ((b + 1) <= (e - 1))
                border[b+1:e-1] = fn(border[b:e])
            end
        end
        return border
    end

    function remove_line(coef::AbstractVecOrMat)
        # 1. in function bdry_recur, input for process_border data[t:b,mc]
        #    is a vector
        # 2. so the coef here is also a vector, so we only consider the length of Tuple
        #    and n is always 1

        m = length(coef)

        n = m
        coef = coef'

        if (m == 2 & n == 2)
            remove_line = coef[2]' .- [1]' / 1 * (coef[2] - coef[1]) .- coef[1]
        else
            remove_line = coef[2:end-1]' - [1:1:n-2;]' / (n - 1) * (coef[end] - coef[1]) .- coef[1]
        end

        remove_line = remove_line'

        return remove_line
    end

    function addb_line(coef::AbstractVecOrMat)
        (m, n) = size(coef)
        if (n == 1)
            do_transf = 1
            coef = coef'
            n = m
        else
            do_transf = 0
        end
        addb_line = coef[2:end-1]' + [1:1:n-2;]' / (n - 1) * (coef[end] - coef[1]) .+ coef[1]
        if (do_transf == 1)
            addb_line = addb_line'
        end
        return addb_line
    end

        
    function gridBlkBdryTrans(data::AbstractArray, ll::AbstractArray, func::Function, inverse::Bool)
        if (inverse === nothing)
            inverse = false
        end

        (m, n) = size(data)
        # Initialize boundary endpoint indices
        #          b4
        #    +------------+
        #    |            |
        #    |            |
        # b1 |     im     | b3
        #    |            |
        #    |            |
        #    +------------+
        #          b2
        b1 = [1 m]
        b2 = [1 n]
        b3 = [1 m]
        b4 = [1 n]
        # Now we find the boundaries and process them.
        ix = 1  # reset index to process borders

        # nested function for boundary recurrent
        function bdry_recur(l::Int, r::Int, t::Int, b::Int, lb1::AbstractArray, lb2::AbstractArray, lb3::AbstractArray, lb4::AbstractArray, level::Int, inverse::Bool)
            if (ll[ix] != level)
                mc = floor(Int, (l + r) / 2)
                mr = floor(Int, (b + t) / 2)

                # Update boundary vectors with new partition
                lb1 = sort(union(lb1, mr))'
                lb2 = sort(union(lb2, mc))'
                lb3 = sort(union(lb3, mr))'
                lb4 = sort(union(lb4, mc))'

                # We get two new boundary vectors from new interior boundaries
                lb6 = [l, mc, r]
                lb5 = [t, mr, b]

                (lb1, lb6, lb5, lb4) = bdry_recur(l, mc, t, mr, lb1, lb6, lb5, lb4, level + 1, inverse)
                (lb1, lb2, lb5, lb6) = bdry_recur(l, mc, mr, b, lb1, lb2, lb5, lb6, level + 1, inverse)
                (lb5, lb6, lb3, lb4) = bdry_recur(mc, r, t, mr, lb5, lb6, lb3, lb4, level + 1, inverse)
                (lb5, lb2, lb3, lb6) = bdry_recur(mc, r, mr, b, lb5, lb2, lb3, lb6, level + 1, inverse)

                lb5 = lb5 .- t .+ 1
                lb6 = lb6 .- l .+ 1

                # data_lb5 = func(data[t:b,mc], lb5, inverse ? ibdrytrans : fbdrytrans);
                # count = 1;
                # for i in t:b
                #     data[1, mc] = data_lb5[count];
                #     count = count + 1;
                # end

                data[t:b, mc] = func(data[t:b, mc], lb5, inverse ? ibdrytrans : fbdrytrans)
                data[mr, l:r] = func(data[mr, l:r], lb6, inverse ? ibdrytrans : fbdrytrans)
            else
                ix = ix + 1
            end

            return (lb1, lb2, lb3, lb4)
        end
        # end nested function

        # start the boundary recurrent for (b1,b2,b3,b4)
        (b1, b2, b3, b4) = bdry_recur(1, n, 1, m, b1, b2, b3, b4, 0, inverse)

        data[:, 1] = func(data[:, 1], b1, inverse ? ibdrytrans : fbdrytrans)
        data[1, :] = func(data[1, :], b4, inverse ? ibdrytrans : fbdrytrans)
        data[:, end] = func(data[:, end], b3, inverse ? ibdrytrans : fbdrytrans)
        data[end, :] = func(data[end, :], b2, inverse ? ibdrytrans : fbdrytrans)
        return data
    end

    function gridBlkTrans(data::AbstractArray, ll::AbstractArray, func::Function)

        # Recursive 2D Block Transform Function


        data = blk_recursive(data, 0, ll, func)

        return data
    end

    function blk_recursive(data::AbstractArray, level::Int, ll::AbstractArray, func::Function)
        ix = 1 # initialize level list index
        if level == ll[ix]
            # note that this only transforms the interior of im. Boundary
            # is left untouched.
            trans = func(data)
            ix = ix + 1
        else
            (mm, nn) = size(data)
            m2 = floor(Int, (mm + 1) / 2) # changed from (mm+1)/2
            n2 = floor(Int, (nn + 1) / 2) # changed from (nn+1)/2
            t1 = blk_recursive(data[1:m2, 1:n2], level + 1, ll, func)
            t2 = blk_recursive(data[m2:end, 1:n2], level + 1, ll, func)
            t3 = blk_recursive(data[1:m2, n2:end], level + 1, ll, func)
            t4 = blk_recursive(data[m2:end, n2:end], level + 1, ll, func)

            # note that the regions share some boundaries, hence we have to
            # be alittle carefull when combining
            trans = [t1[1:end-1, 1:end-1] t3[1:end-1, :]
                t2[:, 1:end-1] t4]
        end
        return trans
    end

end
