# phlct2d.m:
#   A function to calculate the 2d-phlct coefficints of an input 2d-image.
#      (c) 2006 by Katsu Yamatani
#         All rights reserved.
#
# Usage:    V = phlct2d( F, N )
#
#
# Inputs:
#    F  :  m x n matrix.
#          N x N block based 2d-DCT coefficients,
#          Note that m and n must be multiples of N.
#
#    N  :  brock size.
#
# Output:
#    V  :  2d-phlct coefficints V := F-U.
#          m x n matrix.
#

module PHLCT2D

    export phlct2d

    
    function phlct2d( F, N )

        # set eta...
        eta1d, eta2d = set_eta2d( N );

        # compute U component.
        U = SetU1( F, eta1d, N ) + SetU2( F, eta2d, N );

        # 2d-phlct coefficints.
        V = F - U;
        return V;
    end

    # SetU1.m:
    #   A function to calculate the 2d-DCT coefficints of the u component.
    #    Used for `` (k, 0), (0, k), k=1,2,...N-1'' parts.
    #   See eqs. (13) in " Improvement of DCT-based compression algorithms 
    #      using Poisson's equation (Katsu Yamatani and Naoki Saito), 
    #       IEEE Transactions on Image Processing."
    #      (c) 2006 by Katsu Yamatani
    #         All rights reserved.
    #
    # Usage:
    #    U = SetU1( F, eta, N )
    #
    # Inputs:
    #    F  :  m x n matrix.
    #          N x N block based 2d-DCT coefficients,
    #          the number of blocks is "km x kn = m/N x n/N".
    #          Note that m and n must be multiples of N.
    #
    #    N  :  brock size.
    #   eta :  coefficints used in eqs.(13); 
    #            see Table I ( an example with N=8 ).
    #
    #
    # Output:
    #    U  :  2d-DCT coefficints of the u component.
    #     `` (k, 0), (0, k), k=1,2,...N-1'' parts only.
    #
    function SetU1( F, eta, N )

        # Initial parameter and array setup.
        m, n = size( F );
        U = zeros(m,n);
        km = Int(m/N);     # the number of blocks along the x-axis
        kn = Int(n/N);     # the number of blocks along the y-axis
        eta1 = eta[:,1]; #eta1[:,1] = eta;
        eta2 = eta1;
        eta2[1:2:N-1,1] =  eta1[1:2:N-1,1];
        eta2[2:2:N-1,1] = -eta1[2:2:N-1,1]; # reverse

        # Set U on each blocks.
        #   Approximate the normal derivative at block boundaries
        #   by using the mean value of each rows or columns.
        #   In fact, just compute the difference of DCT coefficints
        #   for adjacent blocks.

        # x-axis
        g = zeros(km+1,kn);
        g[2:km,:] = F[N+1:N:end,1:N:end] - F[1:N:end-N,1:N:end];
        etaf = repeat( eta1, 1, kn );
        etab = repeat( eta2, 1, kn );
        for i = 1:km
            m1 = N*i;       # end point of a current brock
            m0 = m1-N+2;    # start point of a current brock
            U[m0:m1,1:N:end] =  U[m0:m1,1:N:end] + etaf.*repeat( g[i,:], 1, N-1 )' .+ etab.*repeat( g[i+1,:], 1, N-1 )';
        end

        # y-axis
        g = zeros(km,kn+1);
        g[:,2:kn] = F[1:N:end,N+1:N:end] - F[1:N:end,1:N:end-N];
        etaf = repeat( eta1',km,1 );
        etab = repeat( eta2',km,1 );
        for j = 1:kn
            n1 = N*j;       # end point of a current brock
            n0 = n1-N+2;    # start point of a current brock
            U[1:N:end,n0:n1] = U[1:N:end,n0:n1] + etaf.*repeat( g[:,j], 1, N-1 ) .+  etab.*repeat( g[:,j+1], 1, N-1 );
        end

        return U;

    end

    # SetU2.m:
    #   A function to calculate the 3d-DCT coefficints of the u component.
    #    Used for `` (i, j), i,j=1,2,...,N-1'' parts.
    #   See eqs. (12) in " Improvement of DCT-based compression algorithms 
    #      using Poisson's equation (Katsu Yamatani and Naoki Saito), 
    #       IEEE Transactions on Image Processing."
    #      (c) 2006 by Katsu Yamatani
    #         All rights reserved.
    #
    # Usage:
    #    U = SetU2( F, eta, N )
    #
    # Inputs:
    #    F  :  m x n matrix.
    #          N x N block based 2d-DCT coefficients,
    #          the number of blocks is "km x kn = m/N x n/N".
    #          Note that m and n must be multiples of N.
    #
    #    N  :  brock size.
    #   eta :  coefficints used in eqs.(12); 
    #            see Table I ( an example with N=8 ).
    #
    #
    # Output:
    #    U  :  2d-DCT coefficints of the u component.
    #     `` (i, j), i,j=1,2,...,N-1'' parts only.
    #
    function SetU2( F, eta, N )

        # Initial parameter and array setup.
        m, n = size( F );
        U = zeros(m,n);
        km = Int(m/N);     # the number of blocks along the x-axis
        kn = Int(n/N);     # the number of blocks along the y-axis
        eta1 = zeros(N-1,N);
        eta1[:,2:N] = eta;
        eta2 = eta1;
        eta2[1:2:N-1,:] =  eta1[1:2:N-1,:]; ##necessary?
        eta2[2:2:N-1,:] = -eta1[2:2:N-1,:]; # reverse
    
        # Set U on each blocks.
        #   Approximate the normal derivative at block boundaries
        #   by using the mean value of each rows or columns.
        #   In fact, just compute the difference of DCT coefficints
        #   for adjacent blocks.
    
    
        # x-axis
        g = zeros(km+1,n);
        g[2:km,:,:] = F[N+1:N:end,:] - F[1:N:end-N,:];
        etaf = repeat( eta1, 1, kn );
        etab = repeat( eta2, 1, kn );
    
        for i = 1:km
        m1 = N*i;       # end point of a current brock
        m0 = m1-N+2;    # start point of a current brock
        U[m0:m1,:] = U[m0:m1,:] + etaf.*repeat( g[i,:], 1, N-1 )' .+ etab.*repeat( g[i+1,:], 1, N-1 )';
        end
    
    
        # y-axis
        g = zeros(m,kn+1);
        g[:,2:kn] = F[:,N+1:N:end] - F[:,1:N:end-N];
        etaf = repeat( eta1',km,1 );
        etab = repeat( eta2',km,1 );
        for j = 1:kn
        n1 = N*j;       # end point of a current brock
        n0 = n1-N+2;    # start point of a current brock
        U[:,n0:n1] = U[:,n0:n1] + etaf.*repeat( g[:,j], 1, N-1 ) .+  etab.*repeat( g[:,j+1], 1, N-1 );
        end
    
        return U;
    end
end