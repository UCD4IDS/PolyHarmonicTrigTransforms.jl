try
    using Pkg
    Pkg.activate(".")
    using Documenter
    @info "Documenter loaded" version=Documenter.VERSION
    include(joinpath(@__DIR__, "..", "docs", "make.jl"))
    @info "make.jl executed"
catch e
    println("ERROR: ", e)
    showerror(stderr, e)
    println()
    for (i, frame) in enumerate(stacktrace(catch_backtrace()))
        println(i, ": ", frame)
    end
    rethrow()
end
