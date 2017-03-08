% SNUFFTD Shifting non-uniform fast Fourier transform (S-NUFFT)
%
% Usage
%    f = snufftd(N, omega, alpha, b, q, m);
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

function f = snufftd(N, omega, alpha, b, q, m, use_mx)
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

    n = size(omega, 1);
    d = size(omega, 2);

    mu = round(m*omega);

    delta = m*omega-mu;

    T = zeros([N*ones(1, d) 1]);

    grid = make_grid(N, d);

    phase_shifts = exp(2*pi*i*grid/(m*N));

    phase_shifts = num2cell(phase_shifts, 1:d);
    phase_shifts = phase_shifts(:);

    % Loop through all combinations of grid shifts, apply sub-NUFFT and
    % aggregate.
    for grid_ind = 1:m^d
        grid_shift = cell(d, 1);
        [grid_shift{:}] = ind2sub(m*ones(1, d), grid_ind);
        grid_shift = cell2mat(grid_shift)-1;

        T_sub = sub_snufftd(N, grid_shift, omega, alpha, b, q, m, use_mx);

        for k = 1:d
            for l = 1:grid_shift(k)
                T_sub = T_sub.*phase_shifts{k};
            end
        end

        T = T + T_sub;
    end

    wt = exp(b*(2*pi*[-N/2:N/2-1]'/(m*N)).^2);

    f = T;
    for l = 1:d
        f = bsxfun(@times, f, permute(wt, [2:l 1 l+1]));
    end
end

function grid = make_grid(N, d)
    rng = arrayfun(@(k)([-k/2:k/2-1]'), N*ones(1, d), 'uniformoutput', false);

    grid = cell(d, 1);
    [grid{:}] = ndgrid(rng{:});

    grid = cell2mat(permute(grid, [2:d+1 1]));
end

function T = sub_snufftd(N, grid_shift, omega, alpha, b, q, m, use_mx)
    d = size(omega, 2);

    if ~use_mx
        tau = sub_snufftd_spread(N, grid_shift, omega, alpha, b, q, m);
    else
        tau = sub_snufftd_spread_mx(N, grid_shift, omega, alpha, b, q, m);
    end

    % Apply IFFT and reweight.
    tau = ifftshift(tau);
    T = N^d*ifftn(tau);
    T = fftshift(T);
end
