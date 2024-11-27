module PolyHarmonicTrigTransforms
    
    using AbstractFFTs, FFTW

    include("dst.jl")
    include("idst.jl")
    include("llst.jl")
    include("illst.jl")
    include("solvelaplace.jl")
    include("llst2d.jl")
    include("phlct.jl")
    include("phlct2d.jl")
 
    using .DST, .IDST, .LLST, .LLST2D, .ILLST, .SOLVELAPLACE, .PHLCT, .PHLCT2D

    export dst, idst, llst, llstapprox2, illst, solvelaplace, phlct, phlct2d
end

#using .PolyHarmonicTrigTransforms