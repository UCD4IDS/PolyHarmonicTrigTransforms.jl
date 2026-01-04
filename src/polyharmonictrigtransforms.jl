# Case-correct package entrypoint for case-sensitive filesystems (CI/macOS)
# This simply includes the existing implementation file to maintain the
# current layout while satisfying `Pkg` expectations for `src/ModuleName.jl`.
include(joinpath(@__DIR__, "PolyHarmonicTrigTransforms.jl"))

# Module entrypoint
module PolyHarmonicTrigTransforms

# Include implementation files (relative to this file's directory)
for fname in (
    "dst.jl", "idst.jl", "llst.jl", "illst.jl", "llst2.jl",
    "llstapprox2.jl", "solvelaplace.jl", "solvelaplace_old.jl",
    "phlct.jl", "phlct2d.jl", "qtllst2dl1.jl", "helper.jl",
)
    include(joinpath(@__DIR__, fname))
end

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
