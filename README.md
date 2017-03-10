# Shifting non-uniform fast Fourier transform in d dimensions (SNUFFTd)

Implements non-uniform fast Fourier transforms following [1, 2] with
modifications to reduce memory usage. Instead of computing one large fast
Fourier transform on the oversampled grid, the SNUFFTd words on shifted
sub-grids. As a result, no arrays larger than the final output have to be
allocated.

### References

[1] A. Dutt and V. Rokhlin, Fast Fourier transforms for nonequispaced data,
SIAM Journal on Scientific Computing, 14 (1993), pp. 1368-1393.

[2] L. Greengard and J.-Y. Lee, Accelerating the nonuniform fast fourier
transform, SIAM review, 46 (2004), pp. 443–454.
