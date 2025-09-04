module PolyHarmonicTrigTransforms
    
    using AbstractFFTs, FFTW

    include("dst.jl")
    include("idst.jl")
    include("llst.jl")
    include("illst.jl")
    include("llst2.jl")
    include("llstapprox2.jl")
    include("solvelaplace.jl")
    include("solvelaplace_old.jl")
    include("phlct.jl")
    include("phlct2d.jl")
    include("qtllst2dl1.jl")

    include("helper.jl")
 
    using 
        .DST,
        .IDST,
        .LLST,
        .ILLST,
        .LLST2,
        .LLSTAPPROX2,
        .SOLVELAPLACE,
        .SOLVELAPLACE_OLD,
        .PHLCT,
        .PHLCT2D,
        .QTLLST2DL1,
        .HELPER

    export 
        #dst.jl
        dst, dst_old,
        idst, idst_old,
        #llst2
        llst2d, illst2d,
        illst,
        llst2,
        llstapprox2,
        solvelaplace,
        solvelaplace_old,
        #phlct,
        phlct_backward, phlct_forward, phlct_restore
        phlct2d,
        qtllst2dl1,
        #helper.jl
        drawpartition2d2, l2norm, is_level_list_valid, is_valid_subband
end

#using .PolyHarmonicTrigTransforms