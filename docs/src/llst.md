llFor Intel might need to run for MKL 
FFTW.set_provider!("mkl")
# DST Discrete Sine Transform.
Computes the type I discrete sine transform of 'x'.  If 'n' is given, 
then 'x' is padded or trimmed to length 'n' before computing the 
transform. If 'x' is a matrix, compute the transform along the columns of 
the matrix.

The discrete sine transform X of x can be defined as follows:
`NX[k] = sum x[n] sin (pi n k / (N+1)),  k = 1, ..., N
    n=1
`
Syntax:
y = DST(x)
y = DST(x, n)

Inputs:
x - vector to be transformed
n - length of tranform to perform

Outputs:  
y - tranformed vector

Author:
Jason Hee Yoon (2022)
Shuchen Ye (2022)

See also: IDST