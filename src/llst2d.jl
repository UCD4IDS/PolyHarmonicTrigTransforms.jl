#
# llstapprox2.m: A function to automatically evaluate the quality of
# LLST-based approximation in terms of PSNR only.
# All the corner points are kept by default.
#
# Usage:
#
# psnr=llstapprox2(data, levlist, krange, bdrykeep)
#
# Inputs: data: Input image array.
#         levlist: levels list array, currently only homogeneous
#                  split works correctly.
#         krange: a vector indices k to be used for nonlinear approximation.
#         bdrykeep: 0 if boundary coefficients are also sorted as a
#                   part of coefficient selection (default).
#                   1 if boundary coefficients are retained and the
#                   only inside coefficients are sorted.
#
#	(c) 2002-2004 PHLTT Development Team @ Math Dept. UC Davis
#			    All rights reserved.
#	                    Coded by Naoki Saito
#

module LLST2D
    include("llst.jl")
    include("illst.jl")
    include("util.jl")
    using .LLST
    using .ILLST
    using .UTIL

    export llstapprox2
    @inline function llstapprox2(data::AbstractArray, levlist::AbstractArray, krange::AbstractArray, bdrykeep::Bool=false)
        # Sanity check.
        (M, N) = size(data)
        round.(data, digits=4)
        #levlist=Float64.(levlist);
        #if (max(krange) > length(data(:))-(M/2^levlist(1)+1)*(N/2^levlist(1)+1))
        #  error(sprintf('max(krange) = #d exceeds the limit!', max(krange)));
        #end
        if (prod(diff(krange)) < 0)
            #error('Elements of krange is not in increasing order.');
        end
        klen = length(krange)
        #levlist=Float64.(levlist);

        # Do the transform first.
        coef = llst(data, levlist)
        (M, N) = size(coef) # just in case, recomputes M, N.

        # Split the coefficients into 'in', 'bo', 'co'.
        (in, bo, co) = split_llst_coefs(coef, levlist)
        if (!isempty(data) && !isempty(levlist) && !isempty(krange) && bdrykeep)
            inbo = in
        else
            a = in
            b = bo

            inbo = [a' b']'
        end
        #disp(sprintf('Number of corner points:         #d\n', length(co(:))));
        #disp(sprintf('Number of boundary coefficients: #d\n', length(bo(:))));
        #disp(sprintf('Number of inside coefficients:   #d\n', length(in(:))));

        # Sort the 'inbo' coefficients.
        ind = sortperm(abs.(Float64.(inbo)))
        #ind = sort(tmp);

        # Prepare the output arrays.
        psnr = zeros(klen, 1)
        mssim = zeros(klen, 1)

        # Start the loop of evaluations.
        fact = sqrt(length(data))
        rmax = maximum(abs.(data))
        if (!isempty(data) && !isempty(levlist) && !isempty(krange) && bdrykeep)
            Threads.@threads for k = 1:klen
                tmp = copy(inbo)
                tmp[ind[1:end-krange(k)]] .= 0
                tmp = merge_llst_coefs(tmp[1:end], bo, co, levlist, M, N)
                #if rem(k,1000)==0
                #disp(sprintf('k=%d',k))
                #end
                tmp = illst(tmp, levlist)
                psnr[k] = 20 * log10(rmax * fact / norm(data - tmp))
            end
        else
            Threads.@threads for k = 1:klen
                @info "k: " k Dates.Time(Dates.now())
                tmp = copy(inbo)
                tmpInd = ind[1:end-krange[k]]
                tmp[tmpInd] .= 0
                tmp = merge_llst_coefs(tmp[1:length(in)], tmp[length(in)+1:end], co, levlist, M, N)
                #if rem(k,1000)==0
                #disp(sprintf('k=%d',k))
                #end
                tm = illst(tmp, levlist)
                tnorm = norm(data - tm, 2)
                if (tnorm <= eps())
                    tnorm = eps()
                end
                psnr[k] = 20 * log10(rmax * fact / tnorm)
            end
        end
        return sort(psnr, dims=1)
    end
end