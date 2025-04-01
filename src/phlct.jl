module PHLCT

    using Statistics, FFTW
    export phlct_backward, phlct_forward, phlct_restore
    
    function phlct_backward(input::AbstractVecOrMat, N::Int)
        # 
        # input
        # input : the m x n DCT coefficients of the residual image
        #  N    : the size of each block
        # 
        # output
        # out: the m x n reconstructed data from DCT coefficients and PHfunction
        # 

        m, n = size(input)
        km = m ÷ N
        kn = n ÷ N

        # reconstruct residual image u from DCT coefficients
        u = zeros(m, n)
        h = pi / (2 * N)
        t = h * (repeat(reshape(collect(0:2*N-1), 1, 2 * N), 2 * N, 1) .+ repeat(reshape(collect(0:2*N-1), 2 * N, 1), 1, 2 * N))
        t = (2 * N)^2 * (cos.(t) .+ im * sin.(t))
        buf = zeros(2 * N, 2 * N)
        z1 = zeros(N, 1)
        z2 = zeros(N - 1, 1)
        z3 = zeros(1, 2 * N)

        for k = 1:km
            m1 = N * k
            m0 = m1 - N + 1
            for i = 1:kn
                n1 = N * i
                n0 = n1 - N + 1
                s = input[m0:m1, n0:n1]
                buf = [s z1 -s[:, N:-1:2]; z3; -s[N:-1:2,:] z2 s[N:-1:2, N:-1:2] ]
                temp = ifft(t .* buf)
                u[m0:m1, n0:n1] = real.(temp[1:N, 1:N])
            end
        end

        # reconstruct PHfunction v from DC components
        dcc = input[1:N:m, 1:N:n]
        v = ndm(dcc, N)

        # add u and v
        out = u + v

        return out
    end

    function phlct_forward(input::AbstractVecOrMat, N::Int)
        # 
        # input
        # input : an m x n image data
        #  N    : the size of each block
        # 
        # output
        # out: m x n DCT coefficients of u := in-v, where v denotes PHfunction
        #  
    
        m, n = size(input)
        km = Int(m / N)
        kn = Int(n / N)
    
        # compute DC components
    
        dcc = zeros(km, kn)
        for k = 1:km
            m1 = N * k
            m0 = m1 - N + 1
            for i = 1:kn
                n1 = N * i
                s = input[m0:m1, n1-N+1:n1]
                dcc[k, i] = mean(s)
            end
        end
    
        # remove PHfunction v from input image data
        input = input - ndm(dcc, N)
    
        # compute DCT coefficients of u = input-v
        h = pi / (2 * N)
        t = h * (repeat(reshape(collect(0:N-1), 1, N), N, 1) .+ repeat(reshape(collect(0:N-1), N, 1), 1, N))
        t = (cos.(t) .- im * sin.(t)) / (2 * N)^2
        buf = zeros(2 * N, 2 * N)
    
        out = zeros(m, n)
        for k = 1:km
            m1 = N * k
            m0 = m1 - N + 1
            for i = 1:kn
                n1 = N * i
                n0 = n1 - N + 1
                s = input[m0:m1, n0:n1]
                buf[1:N, 1:N] = s
                buf[1:N, N+1:end] = s[:, N:-1:1]
                buf[N+1:end, 1:N] = s[N:-1:1, :]
                buf[N+1:end, N+1:end] = s[N:-1:1, N:-1:1]
                temp = fft(buf)
                # temp = fft(buf, [size(buf, 1), size(buf, 2)])
                # display(temp)
                out[m0:m1, n0:n1] .= real.(t .* temp[1:N, 1:N])
            end
        end
        return out
    end
    

    # 
    # phlct_restore.m:
    #     A function to restore the truncated 
    #     DCT coefficints based on PHLCT.
    # 
    # Usage:
    # 
    #     [ out ] = phlct_restore( in, qt )
    # 
    # Inputs
    #    in :  m x n matrix.
    #          N x N block based DCT coefficients
    #          after the JPEG quantization, i.e.,
    #          the number of blocks is "m/N x n/N".
    # 
    #    qt :  N x N matrix.
    #          Quantization table.  
    # 
    # 
    # Output
    #    out:  DCT coefficints which are restored by PHLCT

    function phlct_restore(in::AbstractVecOrMat, qt::AbstractVecOrMat)

        # Set the quantization table for each blocks
        m, n = size(in)
        N = length(qt)
        q = repeat(qt, div(m, N), div(n, N))

        # DCT coefficients of PHfunction
        dcu = dcndm(in, N)

        # Set dcu[i, j] <- 0 if abs(dcu[i, j]) > q[i, j]/2
        dcu = (1 - sign.(abs.(round.(dcu ./ q)))) .* dcu

        # Restoration of truncated coefficients
        # out[i, j] <- dcu[i, j] if in[i, j] == 0
        # out[i, j] <- in[i, j]  else
        out = in + (1 - sign.(abs.(in))) .* dcu

        return out
    end
    
    #########PRIVATE#########

    #
    # ndm.m:
    #   A function to calculate the DCT coefficints (except 0th) of PHLCT function.
    #
    # Usage:
    #
    #     [ out ] = ndm( in, N )
    #
    # Inputs:
    #    in :  m x n matrix.
    #          N x N block based DCT coefficients,
    #          the number of blocks is "m/N x n/N".
    #
    #   Note that m and n must be multiples of N.
    #
    #
    # Output:
    #    out:  DCT coefficints (except 0th) of PHLCT function.
    #
    function ndm( input, N )

        km, kn = size(input)

        out = zeros(N * km, N * kn)
    
        g = zeros(km + 2, kn + 2)
        g[2:km+1, 2:kn+1] = input
        g[1, 2:kn+1] .= input[1, :]
        g[km+2, 2:kn+1] .= input[km, :]
        g[2:km+1, 1] .= input[:, 1]
        g[2:km+1, kn+2] .= input[:, kn]
    
        temp = collect(0:N-1) .- 0.5 * (N - 1)
         x = repeat(temp, 1, N)
         y = x'
        #y = repeat(temp, 1, N)
        #x = y'
    
        for i = 2:km+1
            m1 = N * (i - 1)
            m0 = m1 - N + 1
            for j = 2:kn+1
                n1 = N * (j - 1)
                a = (g[i, j-1] .+ g[i, j+1] .- 2 .* g[i, j]) ./ N
                b = (g[i-1, j] .+ g[i+1, j] .- 2 .* g[i, j]) ./ N
                c = g[i, j+1] .- g[i, j-1]
                d = g[i+1, j] .- g[i-1, j]
                out[m0:m1, n1-N+1:n1] .= ((a .* x.^2 .+ b.* y.^ 2 .+ c .* x .+ d .* y .- (N^2 - 1) .* (a .+ b) ./ 12) ./ (2 * N))
            end
        end
    
        return out
    end

    #
    # dcndm.m:
    #   A function to calculate the DCT coefficints of block based PHLCT function.
    #
    # Usage:
    #
    #     [ out ] = dcndm( in, N )
    #
    # Inputs:
    #    in :  m x n matrix.
    #          N x N block based DCT coefficients,
    #          the number of blocks is "m/N x n/N".
    #
    #   Note that m and n must be multiples of N.
    #
    #
    # Output:
    #    out:  DCT coefficints of block based PHLCT function.
    #
    function dcndm(in, N)

        # Initial parameter and array setup.
        m, n = size(in)
        out = zeros(m, n)
        km = Int(m / N)     # the number of block rows
        kn = Int(n / N)     # the number of block columns
        gv = zeros(m, kn + 1)
        gh = zeros(km + 1, n)
        sv = zeros(N)
        sh = zeros(N)


        # Set PHLCT function (common part to all blocks).
        y = -0.5 .+ collect(1:N)
        Y = repeat(y, 1, N - 1)
        a = repeat((pi / N) .* collect(1:N-1)', N, 1)
        f = (-(y .- N) .^ 2) / (2 * N)    # 0th (quadric function)
        #f[:, 2:N] = -(exp(-a.*Y)+exp(a.*Y-2*N*a))./((1.0-exp(-2*N*a)).*a); # higher-order
        f = -(exp.(-a .* Y) + exp.(a .* Y - 2 * N .* a)) ./ ((1.0 .- exp.(-2 * N .* a)) .* a)

        # DCT coefficints of PHLCT function (common part to all blocks).
        fm = sort(f,dims=1, rev=true)
        fc = fft([f; fm])
        l = (pi / (2 * N)) .* repeat(collect(1:N-1), 1, N)
        ff = real((cos.(l) - 1im .* sin.(l)) .* fc[2:N]) ./ (N^2 * sqrt(2.0))
        fb = ff
        fb[2:2:N-1, :] = -ff[2:2:N-1, :]   # reverse
        w1 = ff    # used for matching at y=0.
        w2 = fb    # used for matching at y=N.
        w3 = ff'   # used for matching at x=0.
        w4 = fb'   # used for matching at x=N.


        # Approximate the normal derivative at block boundaries
        # by using the mean value of each rows or columns.
        # In fact, just compute the difference of 0th DCT coefficints
        # for adjacent blocks.
        gv[:, 2:kn] = in[:, N+1:N:n] - in[:, 1:N:n-N]
        gh[2:km, :] = in[N+1:N:m, :] - in[1:N:m-N, :]


        # Compute DCT coefficints of PHLCT function on each blocks.
        for k = 1:km
            m1 = N * k
            m0 = m1 - N + 1
            for i = 1:kn
                n1 = N * i
                n0 = n1 - N + 1

                # Set boundary conditions.
                g1 = gh[k, n0:n1]
                g2 = gh[k+1, n0:n1]
                g3 = gv[m0:m1, i]
                g4 = gv[m0:m1, i+1]

                # Combine 4 components.
                for p = 1:N-1
                    sv[p+1, :] = w1[p, :] .* g1 .+ w2[p, :] .* g2
                    sh[:, p+1] = w3[:, p] .* g3 .+ w4[:, p] .* g4
                end

                out[m0:m1, n0:n1] = sv .+ sh

            end
        end

        return out
    end



end
