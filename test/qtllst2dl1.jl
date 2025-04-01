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

# qtllst2dl1 TEST ###################################################
# Load the file into a Dict
cd(joinpath(dirname(@__FILE__)))
cd("../data")
matfile = matread("bfo.mat")
# Extract the matrix "bfo.mat" from the matfile
bfo = matfile["bfo"]
# display(bfo)

(liste, cost) = qtllst2dl1(Int64.(bfo), 3)

@info liste
@info cost
endTime = Dates.Time(Dates.now())
@info "end: " endTime - startTime
#####################################################################