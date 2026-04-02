begin
    using Plots
    using Dates
    using Profile

    # Standard libraries

end

include("../src/PolyHarmonicTrigTransforms.jl")

using .PolyHarmonicTrigTransforms

startTime = Dates.Time(Dates.now())
@info "startTime: " startTime


function L2error(initial, backward)
    # initial & backward have the same size
    (m, n) = size(initial)
    num = m * n
    total_error = 0
    for i in 1:m
        for j in 1:n
            total_error += abs(initial[m, n] - backward[m, n]) / abs(initial[m, n])
        end
    end

    error = total_error / num
    return error
end

```
llst test
```
# Summary:
#     9x9 matrix & level list [1 1 1 1] - PASSED
#     9x9 matrix & level list [2 2 2 2 1 1 1] - PASSED
#     9x9 matrix & level list [2 3 3 3 3 2 2 1 1 1] - PASSED
#     17x17 matrix & level list [1 1 1 1] - PASSED


#=
mat_9x9 = llst(h, ll)
#@info "mat_9x9" mat_9x9
mat_9x9_inverse = illst(mat_9x9, ll)
@info "mat_9x9_inverse" mat_9x9_inverse
@info "is same" number(mat_9x9_inverse) == h
mat_9x9_err = L2error(h, mat_9x9_inverse)
@info "mat_9x9_err: " mat_9x9_err
=#
n = 9
x = LinRange(-1, 1, n)
y = LinRange(-1, 1, n)
Gaussian = zeros((length(x), length(y)))

Threads.@threads for i in 1:n
    for j in 1:n
        t = exp(-(x[i] + 1 / 3)^2 - (y[j] + 1 / 3)^2)
        Gaussian[i, j] = t
    end
end

# # test for J = 0
# levlist0 = [0]
# krange0 = [1:50:n^2;]
# psnr0 = llstapprox2(Gaussian, levlist0, krange0)

# test for J = 1
levlist1 = [1 1 1 1 1 1]
krange1 = [1:n^2;]
retain1 = krange1 ./ n^2

#data = llst(Gaussian, levlist1) #passed
#=
=#
psnr = @profile llstapprox2(Gaussian, levlist1, krange1)
@info "psnr: " psnr
p = plot(retain1, psnr, show=true);
endTime = Dates.Time(Dates.now())
Profile.print()
savefig(p, "psnr.png")
@info "end: " endTime
@info "took: " Dates.Time(endTime - startTime);
display(p)

# # test for J = 2
# levlist2 = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2]

#=
mat_9x9_2 = llst(h, ll2)
@info "mat_9x9_2" mat_9x9_2
mat_9x9_2_inverse = illst(mat_9x9_2, ll)
@info "mat_9x9_2_inverse" mat_9x9_2_inverse
mat_9x9_2_err = L2error(h, mat_9x9_2_inverse)
@info "mat_9x9_2_err: " mat_9x9_2_err

mat_9x9_3 = llst(h, ll3)
@info "mat_9x9_3" mat_9x9_3
mat_9x9_3_inverse = illst(mat_9x9_3, ll3)
@info "mat_9x9_3_inverse" mat_9x9_3_inverse
mat_9x9_3_err = L2error(h, mat_9x9_3_inverse)
@info "mat_9x9_3_err: " mat_9x9_3_err
=#
# mat_17x17 = llst(i, ll)
# vscodedisplay(mat_17x17)

#####################################################################