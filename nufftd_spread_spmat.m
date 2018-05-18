% NUFFTD_SPREAD_SPMAT Non-uniform fast Fourier transform spreading matrix
%
% Usage
%    tau_mat = nufftd_spread_spmat(N, omega, b, q, m);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%       dimension is d, this results in N^d entries for f.
%    omega: An array of size d-by-n containing the frequencies at which to
%       compute the transform, where n is the number of nodes and d is the
%       dimension. Each entry must be in the range [-N/2, N/2].
%    b, q, m: The parameters of the NUFFT (default b is 1.5629, q is 28, and
%       m is 2).
%
% Output
%    tau_mat: The spreading matrix taking coefficients alpha and spreading
%       them at the proper locations of the oversampled Fourier transform.

function tau_mat = nufftd_spread_spmat(N, omega, b, q, m)
    n = size(omega, 2);
    d = size(omega, 1);

    % Precalculate kernel weights.
    mu = round(m*omega);

    delta = m*omega-mu;

    P1 = exp(-[-q/2:q/2]'.^2/(4*b));
    P2 = exp(-delta.^2/(4*b));
    P3 = exp(2*delta/(4*b));

    tau_mat = sparse([], [], [], (m*N)^d, n, (q+1)^d*n);

    for k = 1:n
        Pfac = bsxfun(@times, P1, P2(1:d,k)');
        Pfac = Pfac .* bsxfun(@power, P3(1:d,k)', [-q/2:q/2]');

        Pk = Pfac(:,1);
        for l = 2:d
            Pk = bsxfun(@times, Pk, permute(Pfac(:,l), [2:l 1 l+1]));
        end
        Pk = Pk * 1/(2*sqrt(b*pi))^d;
        Pk = Pk(:);

        mu_k = bsxfun(@plus, mu(:,k)', [-q/2:q/2]');
        mu_k = mod(mu_k+m*N/2, m*N) + 1;

        mu_k = num2cell(mu_k, 1);
        grid = cell(1, d);
        [grid{:}] = ndgrid(mu_k{:});
        grid = cat(d+1, grid{:});
        grid = reshape(grid, [(q+1)^d d]);

        tau_ind = (grid-1)*((m*N).^[0:d-1]')+1;

        tau_mat(tau_ind,k) = Pk;
    end
end
