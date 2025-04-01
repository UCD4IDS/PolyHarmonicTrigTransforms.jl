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
  using LinearAlgebra

  export drawpartition2d2, l2norm

  function drawpartition2d2(signal::AbstractVecOrMat, liste::AbstractArray, width=0.8)
    (m, n) = size(signal)

    p = plot([1, n, n, 1, 1], [1, 1, m, m, 1], lineWidth=width, aspect_ratio=:equal, linecolor=:red, legend=false)

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
    
    plot!(p, [pn + n / 2, pn + n / 2], [pm, pm + m], lineWidth=width, linecolor=:red, legend=false)
    plot!(p, [pn, pn + n], [pm + m / 2, pm + m / 2], lineWidth=width, linecolor=:red, legend=false)

    
    newpos = recurspartition(p, pm, pn, m / 2, n / 2, liste, pos, level + 1, width)
    newpos = recurspartition(p, pm + m / 2, pn, m / 2, n / 2, liste, newpos, level + 1, width)
    newpos = recurspartition(p, pm, pn + n / 2, m / 2, n / 2, liste, newpos, level + 1, width)
    newpos = recurspartition(p, pm + m / 2, pn + n / 2, m / 2, n / 2, liste, newpos, level + 1, width)
  end

  function l2norm(original::AbstractVecOrMat, mutated::AbstractVecOrMat) 
    return norm(mutated .- original)/ norm(original)
  end
end