begin
    using Dates
    using Profile
    using Plots
    using MAT
end

include("../src/PolyHarmonicTrigTransforms.jl")

using .PolyHarmonicTrigTransforms

startTime = Dates.Time(Dates.now())
@info "startTime: " startTime



matfile = matread("../data/bfo.mat")
bfo = convert(Matrix{Float64}, matfile["bfo"])
bfo128 = convert(Matrix{Float64}, matfile["bfo128"])
qlum = [16 11 10 16 24 40 51 61
        12 12 14 19 26 58 60 55
        14 13 16 24 40 57 69 56
        14 17 22 29 51 87 80 62
        18 22 37 56 68 109 103 77
        24 35 55 64 81 104 113 92
        49 64 78 87 103 121 120 101
        72 92 95 98 112 100 103 99]
n = 8
#####################################################################
# phlct_forward katsu TEST 
#saito = phlct_forward(bfo128, n)
#katsu = phlct2d(bfo128, n)
#iphlct = iphlct2d(katsu, n)
#restore = phlct_restore(forward, qlum)

#backward = phlct_backward(bfo128, n)
#=
@info "phlct_forward: " forward
@info "phlct_restore: " iphlct
#@info "phlct_backward: " backward
print("phlct_forward TEST")
=#
#####################################################################
# phlct_forward Katsu TEST 
#forward = phlct_forward(bfo128, n)
#=
=#
forward = PHLCT.phlct_forward(bfo128, n)
backward = phlct_backward(forward, n)


@info "phlct_forward: " forward
@info "phlct_backward: " backward
print("SAITO TEST")
#####################################################################
# ndm TEST 
#=
N = 3
cd(joinpath(dirname(@__FILE__)))
cd("..")
matfile = matread("dcc.mat")
dcc = convert(Matrix{Float64}, matfile["dcc"])
out = ndm(dcc, N)
@info "ndm: " out
print("ndm TEST")
=#
# # TEST PASSED

#####################################################################
# phlct_backward TEST 
#=
backward = phlct_backward(bfo, 3)
@info "phlct_backward: " backward
print("phlct_backward TEST")

@info "diff" backward - forward
=#