% NUFFTD_SPMAT Non-uniform fast Fourier transform (NUFFT) with sparse matrices
%
% Usage
%    f = nufftd_spmat(N, omega, alpha, b, q, m, use_mx);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%      dimension is d, this results in N^d entries for f.
%    omega: An array of size d-by-n containing the frequencies at which to
%      compute the transform, where n is the number of nodes and d is the
%      dimension. Each entry must be in the range [-N/2, N/2].
%    alpha: An array of length n containing the coefficients.
%    b, q, m: The parameters of the NUFFT (default b is 1.5629, q is 28, and
%       m is 2).
%    use_mx: Specifies whether to use MEX implementation of the spreading
%       function (default false).
%
% Output
%    f: The NUFFT of the frequencies omega with coefficients alpha. In other
%       words,
%
%          f(n1, ..., nd) = sum_{k=1}^n alpha(k)
%              exp(2*pi*i (n1*omega(1,k) + ... + nd*omega(d,k)/N) )
%
%       to some approximation accuracy determined by b, q, and m.

function f = nufftd_spmat(N, omega, alpha, b, q, m)
    if nargin < 4 || isempty(b)
        b = 1.5629;
    end

    if nargin < 5 || isempty(q)
        q = 28;
    end

    if nargin < 6 || isempty(m)
        m = 2;
    end

    % Necessary, for now.
    if mod(N, 2) ~= 0
        error('N must be even.');
    end

    d = size(omega, 1);

    L = size(alpha, 2);

    % Generate sparse spreading matrix
    tau_mat = nufftd_spread_spmat(N, omega, b, q, m);

    % Apply to alphas
    tau = tau_mat*alpha;

    % Convert to full and reshape
    tau = full(tau);
    tau = reshape(tau, [(m*N)*ones(1, d) L]);

    % Apply FFT, flip, and extract central interval.
    for k = 1:d
        tau = ifftshift(tau, k);
    end

    sub1 = [1 m*N:-1:m*N-N/2+2 N/2+1:-1:2];

    idx.type = '()';

    T = tau;
    for k = 1:d
        T = fft(T, [], k);

        idx.subs = repmat({':'}, 1, d+1);
        idx.subs{k} = sub1;

        T = subsref(T, idx);
    end

    for k = 1:d
        T = fftshift(T, k);
    end

    % Reweight to deconvolve.
    f = T;

    wt = exp(b*(2*pi*[-N/2:N/2-1]'/(m*N)).^2);

    for k = 1:d
        f = bsxfun(@times, f, permute(wt, [2:k 1 k+1]));
    end
end
