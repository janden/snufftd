% SNUFFTD Shifting non-uniform fast Fourier transform (S-NUFFT)
%
% Usage
%    f = snufftd(N, omega, alpha, b, q, m, use_mx);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%       dimension is d, this results in N^d entries for f.
%    omega: An array of size d-by-n containing the frequencies at which to
%       compute the transform, where n is the number of nodes and d is the
%       dimension. Each entry must be in the range [-N/2, N/2].
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

    n = size(omega, 2);
    d = size(omega, 1);

    mu = round(m*omega);

    delta = m*omega-mu;

    T = zeros([N*ones(1, d) 1]);

    grid = make_grid(N, d);

    phase_shifts = exp(2*pi*i*grid/(m*N));

    phase_shifts = num2cell(phase_shifts, 1:d);
    phase_shifts = phase_shifts(:);

    % Precompute certain coefficients needed in sub_snufftd.
    precomp = zeros(2, d, n);

    delta = m*omega-round(m*omega);

    precomp(1,:,:) = exp(-delta.^2/(4*b));
    precomp(2,:,:) = exp(-2*delta/(4*b));

    % Loop through all combinations of grid shifts, apply sub-NUFFT and
    % aggregate.
    for grid_ind = 1:m^d
        grid_shift = cell(d, 1);
        [grid_shift{:}] = ind2sub(m*ones(1, d), grid_ind);
        grid_shift = cell2mat(grid_shift)-1;

        T_sub = sub_snufftd(N, grid_shift, omega, alpha, b, q, m, use_mx, ...
            precomp);

        for k = 1:d
            for l = 1:grid_shift(k)
                T_sub = T_sub.*phase_shifts{k};
            end
        end

        T = T + T_sub;
    end

    % Since we're using FFTs, flip to get the IFFT.
    idx.type = '()';

    for l = 1:d
        idx.subs = repmat({':'}, 1, d);
        idx.subs{l} = [1 N:-1:2];

        T = subsref(T, idx);
    end

    % Reweight to deconvolve.
    wt = exp(b*(2*pi*[-N/2:N/2-1]'/(m*N)).^2);

    f = T;
    for l = 1:d
        f = bsxfun(@times, f, permute(wt, [2:l 1 l+1]));
    end
end

function grid = make_grid(N, d)
    % NOTE: The strange structure of the grid is because we are using FFTs
    % instead of IFFTs (since FFTs are faster) and so we need to take this
    % into account when applying the phase shifts.

    rng = arrayfun(@(k)([-k/2 k/2-1:-1:-k/2+1]'), N*ones(1, d), 'uniformoutput', false);

    grid = cell(d, 1);
    [grid{:}] = ndgrid(rng{:});

    grid = cell2mat(permute(grid, [2:d+1 1]));
end

function T = sub_snufftd(N, grid_shift, omega, alpha, b, q, m, use_mx, precomp)
    n = size(omega, 2);
    d = size(omega, 1);

    if ~use_mx
        tau = sub_snufftd_spread(N, grid_shift, omega, alpha, b, q, m, precomp);
    else
        tau = sub_snufftd_spread_mx(N, grid_shift, omega, alpha, b, q, m, precomp);
    end

    % Apply FFT.
    tau = ifftshift(tau);

    T = fftn(tau);

    T = fftshift(T);
end
