begin
    using Dates
    using Profile
    using Plots
end

include("../src/PolyHarmonicTrigTransforms.jl")

using .PolyHarmonicTrigTransforms
import .PolyHarmonicTrigTransforms.HELPER:drawpartition2d2
startTime = Dates.Time(Dates.now())
@info "startTime: " startTime

#=
cd(joinpath(dirname(@__FILE__)))
matfile = matread("../bfo.mat")
bfo = matfile["bfo"]

p=drawpartition2d2(bfo, bfolist)
display(p)
=#
h = [1 2 3 4
    1 2 3 4
    1 2 3 4
    1 2 3 4] 

h1 =  [3 4 5 1
8 2 5 1
1 9 5 1
3 4 5 1
8 2 5 1
1 9 5 1]

#p = drawpartition2d2(h1, [3 3 3])
p = drawpartition2d2(h, [3 3 3])
savefig(p, "drawpartition2d2.png")

endTime = Dates.Time(Dates.now())
@info "endTime: " endTime
@info "Duration: " endTime - startTime