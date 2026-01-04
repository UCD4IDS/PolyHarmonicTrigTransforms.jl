include("../src/polyharmonictrigtransforms.jl")

using Test, .PolyHarmonicTrigTransforms

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

# DST/IDST round-trip
x = randn(16)
y = dst(x)
x_rec = idst(y)
@test isapprox(x, x_rec; atol=1e-10, rtol=1e-8)

# PHLCT forward/backward round-trip
n = 8
img = randn(n, n)
coeffs = phlct_forward(img, n)
img_rec = phlct_backward(coeffs, n)
@test isapprox(img, img_rec; atol=1e-10, rtol=1e-8)

# In-place `_!` variant tests
# DST!/IDST!
x2 = randn(16)
xp = copy(x2)
yp = dst(copy(x2))
dst!(xp)
@test isapprox(xp, yp; atol=1e-10, rtol=1e-8)

id_in = dst(randn(16))
idp = copy(id_in)
idst!(idp)
@test isapprox(idst(id_in), idp; atol=1e-10, rtol=1e-8)

# PHLCT `_!` variants
img2 = randn(n, n)
coeffs_expected = phlct_forward(img2, n)
img2p = copy(img2)
phlct_forward!(img2p, n)
@test isapprox(img2p, coeffs_expected; atol=1e-10, rtol=1e-8)

coeffs2 = phlct_forward(randn(n, n), n)
coeffs2p = copy(coeffs2)
phlct_backward!(coeffs2p, n)
@test isapprox(phlct_backward(coeffs2, n), coeffs2p; atol=1e-10, rtol=1e-8)

# LLST 2D `_!` variants
data = randn(12, 12)
res_ll = llst2d(copy(data))
data_p = copy(data)
llst2d!(data_p)
@test isapprox(data_p, res_ll; atol=1e-10, rtol=1e-8)

res_ill = illst2d(copy(data))
data_p2 = copy(data)
illst2d!(data_p2)
@test isapprox(data_p2, res_ill; atol=1e-10, rtol=1e-8)

# Basic solvelaplace smoke test if available
if isdefined(PolyHarmonicTrigTransforms, :solvelaplace)
	im = zeros(6, 6)
	im[1, :] .= 1.0
	im[end, :] .= 1.0
	res = solvelaplace(im)
	@test size(res) == size(im)
end

println("All tests passed.")