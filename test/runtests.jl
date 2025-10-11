include("../src/PolyHarmonicTrigTransforms")

using .PolyHarmonicTrigTransforms

#check for dst
has_dst = isdefined(PolyHarmonicTrigTransforms, :dst);
@test has_dst = true;

has_idst = isdefined(PolyHarmonicTrigTransforms, :idst);
@test has_idst = true;

has_phlct_backward = isdefined(PolyHarmonicTrigTransforms, :phlct_backward);
@test has_phlct_backward = true;

has_phlct_forward = isdefined(PolyHarmonicTrigTransforms, :phlct_forward);
@test has_phlct_forward = true;

has_phlct_restore = isdefined(PolyHarmonicTrigTransforms, :phlct_restore);
@test has_phlct_restore = true;

has_phlct2d = isdefined(PolyHarmonicTrigTransforms, :phlct2d);
@test has_phlct2d = true;

has_laplace = isdefined(PolyHarmonicTrigTransforms, :laplace);
@test has_laplace::Bool = true;

print("tests passed");
