module QTLLST2DL1

    using LinearAlgebra

    include("dst.jl")
    include("solvelaplace.jl")

    using .DST
    using .SOLVELAPLACE
    
    export qtllst2dl1

    function qtllst2dl1(signal::AbstractVecOrMat, lmax::Int)
        liste, cost = recurs_testdst2d(signal, 0, lmax)
        return liste, cost
    end

    function recurs_testdst2d(signal::AbstractVecOrMat, level::Int, lmax::Int)

        println(size(signal))  # for fixing: print signal size

        (m, n) = size(signal)

        n1 = Int((n - 1) / 2)
        m1 = Int((m - 1) / 2)
        n2 = Int((n + 1) / 2)
        m2 = Int((m + 1) / 2)

        lap = solvelaplace(signal)
        r = signal[2:end-1, 2:end-1] .- lap
        s = sqrt(4 / (m - 1) / (n - 1))
        r1 = dst(dst(r)')'
        r = r1 .* s
        cost = norm(r, 1) #Normalisation since matlab dst is not unitary

        liste = level

        if (level == lmax)
            return liste, cost
        end


        listec1, costc1 = recurs_testdst2d(signal[1:m1+1, 1:n1+1], level + 1, lmax)
        listec2, costc2 = recurs_testdst2d(signal[m1+1:end, 1:n1+1], level + 1, lmax)
        listec3, costc3 = recurs_testdst2d(signal[1:m1+1, n1+1:end], level + 1, lmax)
        listec4, costc4 = recurs_testdst2d(signal[m1+1:end, n1+1:end], level + 1, lmax)
        
        #Deal with the borders
        b1 = signal[:, n2]
        b2 = signal[m2, :]'

        eb1 = b1[2:m2-1] .- (b1[1] .+ collect(1:m2-2) ./ (m2 - 1) .* (b1[m2] - b1[1]))
        eb1 = dst(eb1) .* sqrt(2 / (m2 - 1))

        eb2 = b1[m2+1:m-1] .- (b1[m2] .+ collect(1:m2-2) ./ (m2 - 1) .* (b1[m] - b1[m2]))
        eb2 = dst(eb2) .* sqrt(2 / (m2 - 1))

        eb3 = b2[2:n2-1]' .- (b2[1] .+ collect(1:n2-2)' ./ (n2 - 1) .* (b2[n2] - b2[1]))
        eb3 = dst(eb3) .* sqrt(2 / (n2 - 1))

        eb4 = b2[n2+1:n-1]' .- (b2[n2] .+ collect(1:n2-2)' ./ (n2 - 1) .* (b2[n] - b2[n2]))
        eb4 = dst(eb4) .* sqrt(2 / (n2 - 1))

        costb = norm([eb1[:]; eb2[:]; eb3[:]; eb4[:]], 1)
        
        if (costc1 + costc2 + costc3 + costc4 + costb + abs(b1[m2]) < cost)
            cost = costc1 + costc2 + costc3 + costc4 + costb + abs(b1[m2])
            liste = [listec1 listec2 listec3 listec4]
        end
        return liste, cost
    end
end
