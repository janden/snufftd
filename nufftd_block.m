% NUFFTD_BLOCK Non-uniform fast Fourier transform (NUFFT) of a matrix
%
% Usage
%    f = nufftd_block(N, omega, alpha, b, q, m, use_mx);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%       dimension is d, this results in N^d entries for f.
%    omega: An array of size d-by-n containing the frequencies at which to
%       compute the transform, where n is the number of nodes and d is the
%       dimension. Each entry must be in the range [-N/2, N/2].
%    alpha: A matrix of size n-by-L, where L is the number of vectors to
%       transform.
%    b, q, m: The parameters of the NUFFT (default b is 1.5629, q is 28, and
%       m is 2).
%    use_mx: Specifies whether to use MEX implementation of the spreading
%       function (default false).
%
% Output
%    f: The NUFFT of the frequencies omega with coefficients alpha. In other
%       words,
%
%          f(n1, ..., nd, ell) = sum_{k=1}^n alpha(k, ell)
%              exp(2*pi*i (n1*omega(1,k) + ... + nd*omega(d,k)/N) )
%
%       to some approximation accuracy determined by b, q, and m.

function f = nufftd_block(N, omega, alpha, b, q, m, use_mex)
    if nargin < 7
        use_mex = [];
    end

    d = size(omega, 1);

    L = size(alpha, 2);

    f = zeros([N^d L]);

    for ell = 1:L
        f(:,ell) = reshape(nufftd(N, omega, alpha(:,ell), b, q, m, use_mex), [N^d 1]);
    end

    f = reshape(f, [N*ones(1, d), L]);
end
