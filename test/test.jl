include("../src/PolyHarmonicTrigTransforms.jl")

using .PolyHarmonicTrigTransforms

#check for dst
has_dst = isdefined(PolyHarmonicTrigTransforms, :dst);
@assert has_dst = true;

has_idst = isdefined(PolyHarmonicTrigTransforms, :idst);
@assert has_idst = true;

has_phlct_backward = isdefined(PolyHarmonicTrigTransforms, :phlct_backward);
@assert has_phlct_backward = true;

has_phlct_forward = isdefined(PolyHarmonicTrigTransforms, :phlct_forward);
@assert has_phlct_forward = true;

has_phlct_restore = isdefined(PolyHarmonicTrigTransforms, :phlct_restore);
@assert has_phlct_restore = true;

has_phlct2d = isdefined(PolyHarmonicTrigTransforms, :phlct2d);
@assert has_phlct2d = true;

print("assertions passed");
