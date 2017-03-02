% SUB_SNUFFTD_SPREAD Shifting non-uniform fast Fourier transform spreading
%
% Usage
%    tau = sub_snufftd_spread(N, grid_shift, omega, alpha, b, q, m);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%       dimension is d, this results in N^d entries for f.
%    grid_shift: A vector of length d with integer values between 0 and m-1
%       specifying how much the subsampling grid of the oversampled Fourier
%       transform has been shifted.
%    omega: An array of size n-by-d containing the frequencies at which to
%       compute the transform, where n is the number of nodes and d is the
%       dimension. Each entry must be in the range [-N/2, N/2].
%    alpha: An array of length n containing the coefficients.
%    b, q, m: The parameters of the NUFFT (default b is 1.5629, q is 28, and
%       m is 2).
%
% Output
%    tau: The oversampled Fourier transform tau, sampled at a grid of
%       resolution N shifted by grid_shift. The original Fourier transform
%       is obtained from the non-uniform discrete Fourier transform after
%       convolution with a Gaussian kernel specified by the b, q, and m
%       parameters.

function tau = sub_snufftd_spread(N, grid_shift, omega, alpha, b, q, m)
    n = size(omega, 1);
    d = size(omega, 2);

    mu = round(m*omega);

    delta = m*omega-mu;

    % Precalculate kernel factors and weighting function.
    P1 = exp(-[-q/2:q/2]'.^2/(4*b));
    P2 = exp(-delta.^2/(4*b));
    P3 = exp(2*delta/(4*b));

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
end
