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
%
% Output
%    f: The NUFFT of the frequencies omega with coefficients alpha. In other
%       words,
%
%          f(n1, ..., nd) = sum_{k=1}^n alpha(k)
%              exp(2*pi*i (n1*omega(k,1) + ... + nd*omega(k,d)/N) )
%
%       to some approximation accuracy determined by b, q, and m.

function f = snufftd(N, omega, alpha, b, q, m)
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

    delta = m*omega-mu;

    % Precalculate kernel factors and weighting function. These are sent to
    % each of the sub-NUFFTs.
    P1 = exp(-[-q/2:q/2]'.^2/(4*b));
    P2 = exp(-delta.^2/(4*b));
    P3 = exp(2*delta/(4*b));

    wt = exp(b*(2*pi*[-N/2:N/2-1]'/(m*N)).^2);

    common.mu = mu;
    common.P1 = P1;
    common.P2 = P2;
    common.P3 = P3;
    common.wt = wt;

    f = zeros([N*ones(1, d) 1]);

    grid = make_grid(N, d);

    % Loop through all combinations of grid shifts, apply sub-NUFFT and
    % aggregate.
    for grid_ind = 1:m^d
        grid_shift = cell(d, 1);
        [grid_shift{:}] = ind2sub(m*ones(1, d), grid_ind);
        grid_shift = cell2mat(grid_shift)-1;

        f_sub = sub_snufftd(N, grid_shift, omega, alpha, b, q, m, common);

        phase_shift = reshape(grid, [N^d d])*(grid_shift(:)/(m*N));
        phase_shift = reshape(phase_shift, [N*ones(1, d) 1]);

        f = f + f_sub.*exp(2*pi*i*phase_shift);
    end
end

function grid = make_grid(N, d)
    rng = arrayfun(@(k)([-k/2:k/2-1]'), N*ones(1, d), 'uniformoutput', false);

    grid = cell(d, 1);
    [grid{:}] = ndgrid(rng{:});

    grid = cell2mat(permute(grid, [2:d+1 1]));
end

function f_sub = sub_snufftd(N, grid_shift, omega, alpha, b, q, m, common)
    n = size(omega, 1);
    d = size(omega, 2);

    mu = common.mu;

    delta = m*omega-mu;

    P1 = common.P1;
    P2 = common.P2;
    P3 = common.P3;

    wt = common.wt;

    % Depending on the relative position of each mu frequency on the shifted
    % grid, different ranges of j indices will be needed for the kernel. I.e.,
    % for q = 4 and m = 2, we need [-4 -2 0 2 4] and [-3 -1 1 3]. Enumerate
    % them here.
    js_rngs = arrayfun( ...
        @(shift)([flipud([-shift:-m:-q/2]'); [-shift+m:m:q/2]']), ...
        0:m-1, 'uniformoutput', false);

    % Fill in shifted Fourier transform tau.
    tau = zeros([N*ones(1, d) 1]);

    js_ind = cell(d, 1);

    for k = 1:n
        % Determine the relative shift with respect to the grid and generate
        % the necessary values of j.
        mu_shift = mod(mu(k,:)-grid_shift.', m);

        js_rng = js_rngs(mu_shift+1);

        js_sz = cellfun(@numel, js_rng);

        js_count = prod(js_sz);

        % Look through the j indices and calculate the kernel values at the
        % desired points in tau.
        for j_ind = 1:js_count
            [js_ind{:}] = ind2sub(js_sz, j_ind);
            js = zeros(1, d);
            for l = 1:d
                js(l) = js_rng{l}(js_ind{l});
            end

            mu_j = mod((mu(k,:)+js-grid_shift.')/m+N/2, N)+1;

            tau_ind = (mu_j-1)*(N.^[0:d-1]')+1;

            Pkj = 1/(2*sqrt(b*pi))^d*...
                prod(P1(js+q/2+1).'.*P2(k,1:d).*(P3(k,1:d).^js));

            tau(tau_ind) = tau(tau_ind) + Pkj*alpha(k);
        end
    end

    % Apply IFFT and reweight.
    tau = ifftshift(tau);
    T = N^d*ifftn(tau);
    T = fftshift(T);

    f_sub = T;
    for l = 1:d
        f_sub = bsxfun(@times, f_sub, permute(wt, [2:l 1 l+1]));
    end
end
