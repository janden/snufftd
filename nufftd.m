% NUFFTD Non-uniform fast Fourier transform (NUFFT)
%
% Usage
%    f = nufftd(N, omega, alpha, b, q, m, use_mx);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%      dimension is d, this results in N^d entries for f.
%    omega: An array of size n-by-d containing the frequencies at which to
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
%              exp(2*pi*i (n1*omega(k,1) + ... + nd*omega(k,d)/N) )
%
%       to some approximation accuracy determined by b, q, and m.

function f = nufftd(N, omega, alpha, b, q, m, use_mx)
    if nargin < 4 || isempty(b)
        b = 1.5629;
    end

    if nargin < 5 || isempty(q)
        q = 28;
    end

    if nargin < 6 || isempty(m)
        m = 2;
    end

    if nargin < 7 || isempty(use_mx)
        use_mx = false;
    end

    % Necessary, for now.
    if mod(N, 2) ~= 0
        error('N must be even.');
    end

    d = size(omega, 2);

    % Spread nodes onto oversampled Fourier transform tau.
    if ~use_mx
        tau = nufftd_spread(N, omega, alpha, b, q, m);
    else
        tau = nufftd_spread_mx(N, omega, alpha, b, q, m);
    end

    % Apply IFFT, extract central interval, and reweight.
    tau = ifftshift(tau);
    T = (m*N)^d*ifftn(tau);
    T = fftshift(T);

    sub = {floor(m*N/2)+1-floor(N/2):floor(m*N/2)+ceil(N/2)};
    sub = sub(ones(1, d));
    sub_grid = cell(1, d);
    [sub_grid{:}] = ndgrid(sub{:});
    ind = sub2ind((m*N)*ones(1, d), sub_grid{:});

    f = T(ind);

    wt = exp(b*(2*pi*[-N/2:N/2-1]'/(m*N)).^2);

    for l = 1:d
        f = bsxfun(@times, f, permute(wt, [2:l 1 l+1]));
    end
end
