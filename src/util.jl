module UTIL

    export split_llst_coefs, merge_llst_coefs
    # SPLIT_LLST_COEFS splits the given LLST coefficients expanded with
    # levels list, into a set of vectors representing 'inside',
    # 'border', and 'corner'.  Note that the corner values are pixel
    # values whereas the inside and borders values are DST coefficients
    # (2D, 1D respectively).
    #
    # function [inside,borders,corners]=split_llst_coefs(coef, levlist)
    #
    # Inputs:
    #       coef: a LLST coefficient matrix.
    #       levlist: levels list array used for the DCT expansion.
    # Outputs:
    #        inside: an array for the 2D DST coefficients of the insides.
    #       borders: an array for the 1D DST coefficients of the borders.
    #       corners: an arrary for the corner pixel values.
    #
    #   See also: MERGE_LLST_COEFS
    #
    #	(c) 2002-2004 PHLTT Development Team @ Math Dept. UC Davis
    #		          All rights reserved.
    #                      Coded by Jean-Francois Remy
    #

    @inline function split_llst_coefs(coef::AbstractVecOrMat, levlist::AbstractArray)
        (m, n) = size(coef)

        b1 = [1 n]
        b2 = [1 m]
        b3 = [1 n]
        b4 = [1 m]

        (inside, newindex) = recurs_inside(coef, levlist, 1, 0)
        (borders, corners, newindex, b1, b2, b3, b4) = recurs_borders(coef, 1, n, 1, m, b1, b2, b3, b4, levlist, 1, 0)

        (tmpb, tmpc) = process_border(coef[:, 1], b1)
        borders = [borders; tmpb]
        corners = [corners; tmpc]
        (tmpb, tmpc) = process_border(coef[end, :], b2)
        borders = [borders; tmpb]
        corners = [corners; tmpc[2:end-1]']
        (tmpb, tmpc) = process_border(coef[:, end], b3)
        borders = [borders; tmpb]
        corners = [corners; tmpc]
        (tmpb, tmpc) = process_border(coef[1, :], b4)
        borders = [borders; tmpb]
        corners = [corners; tmpc[2:end-1]']

        return (inside, borders, corners)
    end

    @inline function recurs_inside(coef::AbstractVecOrMat, levellist::AbstractArray, index::Int, level::Int)
        (m, n) = size(coef)
        m2 = trunc(Int, (m + 1) / 2)
        n2 = trunc(Int, (n + 1) / 2)

        if (levellist[index] == level)
            newindex = index + 1
            out = coef[2:end-1, 2:end-1]
            out = out[:]
        else
            (out, newindex) = recurs_inside(coef[1:m2, 1:n2], levellist, index, level + 1)

            (tmp, newindex) = recurs_inside(coef[m2:m, 1:n2], levellist, newindex, level + 1)
            out = [out; tmp]
            (tmp, newindex) = recurs_inside(coef[1:m2, n2:n], levellist, newindex, level + 1)
            out = [out; tmp]
            (tmp, newindex) = recurs_inside(coef[m2:m, n2:n], levellist, newindex, level + 1)
            out = [out; tmp]
        end

        return (out, newindex)
    end

    @inline function recurs_borders(im::AbstractVecOrMat, l::Int, r::Int, t::Int, b::Int, lb1::AbstractVecOrMat, lb2::AbstractVecOrMat, lb3::AbstractVecOrMat, lb4::AbstractVecOrMat, levellist::AbstractArray, index::Int, level::Int)

        (m, n) = size(im)
        m2 = trunc(Int, (m + 1) / 2)
        n2 = trunc(Int, (n + 1) / 2)
        mc = trunc(Int, (l + r) / 2)
        mr = trunc(Int, (b + t) / 2)


        if (levellist[index] != level)
            lb1 = sort(union(lb1, mr))
            lb2 = sort(union(lb2, mc))
            lb3 = sort(union(lb3, mr))
            lb4 = sort(union(lb4, mc))
            lb6 = [l, mc, r]
            lb5 = [t, mr, b]

            (borders, corners, newindex, lb1, lb6, lb5, lb4) = recurs_borders(im[1:m2, 1:n2], l, mc, t, mr, lb1, lb6, lb5, lb4, levellist, index, level + 1)
            (tmpb, tmpc, newindex, lb1, lb2, lb5, lb6) = recurs_borders(im[m2:m, 1:n2], l, mc, mr, b, lb1, lb2, lb5, lb6, levellist, newindex, level + 1)
            borders = [borders; tmpb]
            corners = [corners; tmpc]

            (tmpb, tmpc, newindex, lb5, lb6, lb3, lb4) = recurs_borders(im[1:m2, n2:n], mc, r, t, mr, lb5, lb6, lb3, lb4, levellist, newindex, level + 1)
            borders = [borders; tmpb]
            corners = [corners; tmpc]

            (tmpb, tmpc, newindex, lb5, lb2, lb3, lb6) = recurs_borders(im[m2:m, n2:n], mc, r, mr, b, lb5, lb2, lb3, lb6, levellist, newindex, level + 1)
            borders = [borders; tmpb]
            corners = [corners; tmpc]

            lb5 = lb5 .- t .+ 1
            lb6 = lb6 .- l .+ 1

            (tmpb, tmpc) = process_border(im[:, n2], lb5)
            borders = [borders; tmpb]
            corners = [corners; tmpc[2:end-1]]

            (tmpb, tmpc) = process_border(im[m2, :], lb6)
            ind = findall(x -> x == n2, lb6)[1]
            deleteat!(tmpc, ind) #remove the center corner or else it will appear twice.
            #tmpc[ind] = 0;
            borders = [borders; tmpb]
            if length(tmpc[2:end-1]') > 0
                corners = [corners; tmpc[2:end-1]']
            else
                corners = [corners; []]
            end

            nb1 = lb1
            nb2 = lb2
            nb3 = lb3
            nb4 = lb4
        else
            borders = []
            corners = []
            nb1 = lb1
            nb2 = lb2
            nb3 = lb3
            nb4 = lb4
            newindex = index + 1
        end

        return (borders, corners, newindex, nb1, nb2, nb3, nb4)
    end

    @inline function process_border(border::AbstractVecOrMat, indexes::AbstractArray)
        sort!(indexes)
        n = length(indexes)
        l = length(border)
        corners = border[indexes]
        i = [1:l;]
        #i[indexes] .= 0;
        deleteat!(i, indexes)
        borders = border[i]
        return (borders, corners)
    end

    # MERGE_LLST_COEFS merges the vectors of 'inside', 'border', and
    # 'corner' coefficients for a given levels list into a proper LLST
    # coefficient matrix.
    #
    # function coef = merge_llst_coefs(inside,borders,corners,levlist,m,n)
    #
    # Inputs:
    #        inside: an array for the 2D DST coefficients of the insides.
    #       borders: an array for the 1D DST coefficients of the borders.
    #       corners: an arrary for the corner pixel values.
    #       levlist: a levels list array used for the DCT expansion.
    #             m: a number of rows (fast directions) of the original data.
    #             n: a number of cols (slow directions) of the original data.
    # Outputs:
    #       coef: a LLST coefficient matrix.

    #   See also: SPLIT_LLST_COEFS
    #
    #	(c) 2002-2004 PHLTT Development Team @ Math Dept. UC Davis
    #		          All rights reserved.
    #                      Coded by Jean-Francois Remy
    #
    function merge_llst_coefs(inside::AbstractVecOrMat, borders::AbstractVecOrMat, corners::AbstractVecOrMat, levlist::AbstractArray, m, n)

        b1 = [1 n]
        b2 = [1 m]
        b3 = [1 n]
        b4 = [1 m]

        (coef, newindex, newindinside, newindborder, newindcorner, b1, b2, b3, b4) = 
        recurs(inside, 1, borders, 1, corners, 1, 1, n, 1, m, b1, b2, b3, b4, levlist, 1, 0)

        coef[b1, 1] = corners[newindcorner:newindcorner+length(b1)-1]
        newindcorner = newindcorner + length(b1)
        coef[m, b2[2:end-1]] = corners[newindcorner:newindcorner+length(b2)-2-1]'
        newindcorner = newindcorner + length(b2) - 2
        coef[b3, end] = corners[newindcorner:newindcorner+length(b3)-1]
        newindcorner = newindcorner + length(b3)
        coef[1, b4[2:end-1]] = corners[newindcorner:newindcorner+length(b4)-2-1]'
        newindcorner = newindcorner + length(b4) - 2

        i = [1:m;]
        #i[b1] .= 0;
        b1s = sort(b1)
        toremove = findall(x -> x in b1s, i)
        if (length(toremove) > 0)
            deleteat!(i, toremove)
        end
        coef[i, 1] = borders[newindborder:newindborder+length(i)-1]
        newindborder = newindborder + length(i)
        i = [1:n;]
        #i[b2] .= 0;
        b2s = sort(b2)
        toremove = findall(x -> x in b2s, i)
        if (length(toremove) > 0)
            deleteat!(i, toremove)
        end
        coef[end, i] = borders[newindborder:newindborder+length(i)-1]'
        newindborder = newindborder + length(i)
        i = [1:m;]
        #i[b3] .= 0;
        b3s = sort(b3)
        toremove = findall(x -> x in b3s, i)
        if (length(toremove) > 0)
            deleteat!(i, toremove)
        end
        coef[i, end] = borders[newindborder:newindborder+length(i)-1]
        newindborder = newindborder + length(i)
        i = [1:n;]
        #i[b4] .= 0;
        b4s = sort(b4)
        toremove = findall(x -> x in b4s, i)
        if (length(toremove) > 0)
            deleteat!(i, toremove)
        end
        coef[1, i] = borders[newindborder:newindborder+length(i)-1]'
        newindborder = newindborder + length(i)

        return coef
    end

    #recurs(inside, 1, borders, 
    #1, corners, 1, 
    #1, n, 1, m, 
    #b1, b2, b3, b4, 
    #levlist, 1, 0)
    function recurs(inside::AbstractVecOrMat, indinside::Int, borders::AbstractVecOrMat, 
        indborder::Int, corners::AbstractVecOrMat, indcorner::Int, 
        l::Int, r::Int, t::Int, b::Int, 
        lb1::AbstractVecOrMat, lb2::AbstractVecOrMat, lb3::AbstractVecOrMat, lb4::AbstractVecOrMat, 
        levellist::AbstractArray, index::Int, level::Int)
        m = b - t + 1
        n = r - l + 1
        m2 = trunc(Int, (m + 1) / 2)
        n2 = trunc(Int, (n + 1) / 2)
        mc = trunc(Int, (l + r) / 2)
        mr = trunc(Int, (b + t) / 2)

        out = zeros(m, n)

        if (levellist[index] != level)
            lb1 = sort(union(lb1, mr))
            lb2 = sort(union(lb2, mc))
            lb3 = sort(union(lb3, mr))
            lb4 = sort(union(lb4, mc))
            lb6 = [l, mc, r]
            lb5 = [t, mr, b]

            (out[1:m2, 1:n2], newindex, newindinside, newindborder, newindcorner, lb1, lb6, lb5, lb4) = recurs(inside, indinside, borders, indborder, corners, indcorner, l, mc, t, mr, lb1, lb6, lb5, lb4, levellist, index, level + 1)
            (out[m2:end, 1:n2], newindex, newindinside, newindborder, newindcorner, lb1, lb2, lb5, lb6) = recurs(inside, newindinside, borders, newindborder, corners, newindcorner, l, mc, mr, b, lb1, lb2, lb5, lb6, levellist, newindex, level + 1)
            (out[1:m2, n2:end], newindex, newindinside, newindborder, newindcorner, lb5, lb6, lb3, lb4) = recurs(inside, newindinside, borders, newindborder, corners, newindcorner, mc, r, t, mr, lb5, lb6, lb3, lb4, levellist, newindex, level + 1)
            (out[m2:end, n2:end], newindex, newindinside, newindborder, newindcorner, lb5, lb2, lb3, lb6) = recurs(inside, newindinside, borders, newindborder, corners, newindcorner, mc, r, mr, b, lb5, lb2, lb3, lb6, levellist, newindex, level + 1)

            lb5 = lb5 .- t .+ 1
            lb6 = lb6 .- l .+ 1

            i = [1:m;]
            #i[lb5] .= 0;
            deleteat!(i, lb5)
            out[i, n2] = borders[newindborder:newindborder+length(i)-1]
            newindborder = newindborder + length(i)
            out[lb5[2:end-1], n2] = corners[newindcorner:newindcorner+length(lb5)-2-1]
            newindcorner = newindcorner + length(lb5) - 2

            i = [1:n;]
            #i[lb6] .= 0;
            deleteat!(i, lb6)
            out[m2, i] = borders[newindborder:newindborder+length(i)-1]'
            newindborder = newindborder + length(i)
            ind = findall(x -> x == n2, lb6)[1]
            #lb6[ind] .= 0;
            deleteat!(lb6, ind)
            out[m2, lb6[2:end-1]] = corners[newindcorner:newindcorner+length(lb6)-2-1]'
            newindcorner = newindcorner + length(lb6) - 2


            nb1 = lb1
            nb2 = lb2
            nb3 = lb3
            nb4 = lb4
        else
            b = indinside + (m - 2) * (n - 2) - 1
            a = inside[indinside:b]
            out[2:end-1, 2:end-1] = reshape(a, m - 2, n - 2)
            nb1 = lb1
            nb2 = lb2
            nb3 = lb3
            nb4 = lb4
            newindinside = indinside + (m - 2) * (n - 2)
            newindborder = indborder
            newindcorner = indcorner
            newindex = index + 1
        end

        return (out, newindex, newindinside, newindborder, newindcorner, nb1, nb2, nb3, nb4)
    end
end;