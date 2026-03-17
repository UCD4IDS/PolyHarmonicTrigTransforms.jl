"""
Package entry wrapper to satisfy Pkg's expected source filename.
This file simply includes the actual implementation file `polyharmonictrigtransforms.jl`.
"""
# Module entrypoint
"""
PolyHarmonicTrigTransforms

Provides polyharmonic and trigonometric transform routines:
DST/IDST, LLST families, polyharmonic LCT, and Laplace solvers.

See the `docs/` directory for usage examples and API reference.
"""
module PolyHarmonicTrigTransforms

include("dst.jl")
include("idst.jl")
include("llst.jl")
include("illst.jl")
include("llst2.jl")
include("llstapprox2.jl")
include("phlct.jl")
include("phlct2d.jl")
include("qtllst2dl1.jl")
include("solvelaplace.jl")
include("solvelaplace_old.jl")
include("helper.jl")

using .DST
using .IDST
using .LLST
using .ILLST
using .LLST2
using .LLSTAPPROX2
using .SOLVELAPLACE
using .SOLVELAPLACE_OLD
using .PHLCT
using .PHLCT2D
using .QTLLST2DL1
using .HELPER

# Deprecation wrappers: keep historical `_old` functions but warn and forward
# to the stable APIs. These accept arbitrary args/kwargs to remain compatible.
dst_old(args...; kwargs...) = (Base.depwarn("`dst_old` is deprecated — use `dst`.", :dst_old); dst(args...; kwargs...))
idst_old(args...; kwargs...) = (Base.depwarn("`idst_old` is deprecated — use `idst`.", :idst_old); idst(args...; kwargs...))
solvelaplace_old(args...; kwargs...) = (Base.depwarn("`solvelaplace_old` is deprecated — use `solvelaplace`.", :solvelaplace_old); solvelaplace(args...; kwargs...))

# Explicit public API (keep exports focused and stable)
export dst, dst_old, dst!, plan_dst, idst, idst_old, idst!, plan_idst,
       llst2d, illst2d, llst2d!, illst2d!, illst, llst2, llst2!, llstapprox2,
       solvelaplace, solvelaplace_old,
       phlct_backward, phlct_forward, phlct_backward!, phlct_forward!, phlct_restore,
       phlct2d, qtllst2dl1,
       drawpartition2d2, l2norm, is_level_list_valid, is_valid_subband

end
