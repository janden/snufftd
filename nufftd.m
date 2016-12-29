% NUFFTD Non-uniform fast Fourier transform (NUFFT)
%
% Usage
%    f = nufftd(N, omega, alpha, b, q, m);
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
%
% Output
%    f: The NUFFT of the frequencies omega with coefficients alpha. In other
%       words,
%
%          f(n1, ..., nd) = sum_{k=1}^n alpha(k)
%              exp(2*pi*i (n1*omega(k,1) + ... + nd*omega(k,d)/N) )
%
%       to some approximation accuracy determined by b, q, and m.

function f = nufftd(N, omega, alpha, b, q, m)
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

    n = size(omega, 1);
    d = size(omega, 2);

    mu = round(m*omega);

    % Precalculate kernel factors and weighting function.
    delta = m*omega-mu;

    P1 = exp(-[-q/2:q/2]'.^2/(4*b));
    P2 = exp(-delta.^2/(4*b));
    P3 = exp(2*delta/(4*b));

    wt = exp(b*(2*pi*[-N/2:N/2-1]'/(m*N)).^2);

    % Fill in oversampled Fourier transform tau.
    tau = zeros([m*N*ones(1, d) 1]);

    js_ind = cell(d, 1);

    for k = 1:n
        for j = 1:(q+1)^d
            [js_ind{:}] = ind2sub((q+1)*ones(1, d), j);
            js = zeros(1, d);
            for l = 1:d
                js(l) = js_ind{l}-(q/2+1);
            end

            mu_j = mod(mu(k,:)+js+m*N/2, m*N)+1;

            tau_ind = (mu_j-1)*((m*N).^[0:d-1]')+1;

            Pkj = 1/(2*sqrt(b*pi))^d* ...
                prod(P1(js+q/2+1).'.*P2(k,1:d).*(P3(k,1:d).^js));

            tau(tau_ind) = tau(tau_ind) + Pkj*alpha(k);
        end
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

    for l = 1:d
        f = bsxfun(@times, f, permute(wt, [2:l 1 l+1]));
    end
end
