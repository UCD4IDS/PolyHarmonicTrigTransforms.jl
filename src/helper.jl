# DRAWPARTITION2D2 Display the quadtree decomposition of an image.
#
#   DRAWPARTITION2D2(IM,LEVEL,LINECOL,WIDTH) displays the quadtree 
#   decomposition of IM using the levellist in LEVEL. This will
#   draw the boundaries of each block of the decomposition.
#   The color and width of the boundary lines are controlled by
#   LINECOL ('r' is default) and WIDTH (0.8 is default).
#
#   This function draws the boundary more precisely for the images
#   of size (2^m+1)x(2^n+1), e.g., 129x129 than DRAWPARTITION2D does.
#
#   Note that this function assumes that the image is already drawn
#   before calling this function!
#
#   Remark: This function will not work properly if the dimensions of IM
#   do not allow the quadtree decomposition.
#
#   See also DRAWPARTITION1D, DRAWPARTITION2D.

module HELPER
  
  using Plots
  using FileIO
  using Images
  using LinearAlgebra

  export drawpartition2d2, l2norm, is_level_list_valid, is_valid_subband

  function drawpartition2d2(signal::AbstractMatrix, liste::AbstractMatrix; width=nothing, image=nothing,  fit=false)
    (m, n) = size(signal)

    is_list_valid = is_valid_subband(liste)
    if(!is_list_valid)
      throw(is_list_valid)
    end

    p = plot()
    if(isnothing(width))
      width = 0.8
    end
    if(!isnothing(image))
      if(isfile(image)) 
        img = FileIO.load(image)
        if(fit)
          img = imresize(img, (n, m));
        end
        plot!(p, img)
      else
        print("file not found: ", image)
      end
    end
    plot!(p, [1, n, n, 1, 1], [1, 1, m, m, 1], lineWidth=width, aspect_ratio=:equal, linecolor=:red, legend=false)

    recurspartition(p, 1, 1, m - 1, n - 1, liste, 1, 0, width)

    return p
  end


  function recurspartition(p, pm, pn, m, n, liste::AbstractArray, pos, level::Int, width)
    #@assert length(liste) < pos
    if (isnothing(pos) || pos > length(liste)) 
      return;
    end
    if (level == liste[pos])
      newpos = pos + 1;
      return newpos;
    end
    
    plot!(p, [pn + n / 2, pn + n / 2], [pm, pm + m], linecolor=:red)
    plot!(p, [pn, pn + n], [pm + m / 2, pm + m / 2], linecolor=:red)

    
    newpos = recurspartition(p, pm, pn, m / 2, n / 2, liste, pos, level + 1, width)
    newpos = recurspartition(p, pm, pn + n / 2, m / 2, n / 2, liste, newpos, level + 1, width)
    newpos = recurspartition(p, pm + m / 2, pn, m / 2, n / 2, liste, newpos, level + 1, width)
    newpos = recurspartition(p, pm + m / 2, pn + n / 2, m / 2, n / 2, liste, newpos, level + 1, width)
  end

  function l2norm(original::AbstractVecOrMat, mutated::AbstractVecOrMat) 
    return norm(mutated .- original)/ norm(original)
  end

  function is_level_list_valid(signal::AbstractMatrix, list::AbstractMatrix)
    m,n = size(signal)

    if(!iseven(m) || !iseven(n))
      return false
    end
    return is_valid_subband(list)

  end

  """
  Validate a wavelet-packet levels list encoded as depths of *leaves* in
  depth-first encounter order (UL, UR, LL, LR). Accepts a row vector or vector.

  Rules:
  - The list contains leaf levels only.
  - If a leaf appears deeper than the current expected level, that implies
    recursive subdivision; each subdivision creates 4 children, of which we
    visit the first now and push the other 3 siblings to visit later.
  - The traversal must exactly consume all tokens with no leftovers.

  Returns (Bool, String).
  """
  function is_valid_subband(levels::AbstractVecOrMat{<:Integer})
    a = Int.(vec(levels))
    if isempty(a)
        print("empty input")
        return false
    end

    # Choose a root not deeper than any leaf; picking (min-1) is safe and generic.
    root_level = minimum(a) - 1

    stack = Int[root_level]   # expected subtree roots to visit next (LIFO for DFS)
    i = 1                     # index into a

    while !isempty(stack)
        if i > length(a)
            print("ran out of tokens with $(length(stack)) pending regions")
            return false
        end

        expected = pop!(stack)    # level of the subtree root we are about to visit
        d = a[i]                  # depth of the next leaf token

        if d < expected
            print("token at index $i has level $d < expected $expected")
            return false
        end

        # For each refinement from expected -> d, push the other 3 children
        # (we descend into the first child immediately; siblings are visited later).
        for L in expected:(d-1)
            push!(stack, L+1); push!(stack, L+1); push!(stack, L+1)
        end

        # Consume this leaf
        i += 1
    end

    if i <= length(a)
        print("extra tokens starting at index $i")
        return false
    end

    return true
  end



end