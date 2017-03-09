% NUDFTD Non-uniform discrete Fourier transform (NUFFT)
%
% Usage
%    f = nudftd(N, omega, alpha);
%
% Input
%    N: Desired resolution of the output f, i.e., the side length. If the
%      dimension is d, this results in N^d entries for f.
%    omega: An array of size d-by-n containing the frequencies at which to
%      compute the transform, where n is the number of nodes and d is the
%      dimension. Each entry must be in the range [-N/2, N/2].
%    alpha: An array of length n containing the coefficients.
%
% Output
%    f: The NUDFT of the frequencies omega with coefficients alpha. In other
%       words,
%
%          f(n1, ..., nd) = sum_{k=1}^n alpha(k)
%              exp(2*pi*i (n1*omega(k,1) + ... + nd*omega(k,d)/N) )


function f = nudftd(N, omega, alpha)
    d = size(omega, 1);

    % Necessary, for now.
    if mod(N, 2) ~= 0
        error('N must be even.');
    end

    rng = {[-N/2:N/2-1]'};
    grid = cell(d, 1);
    [grid{:}] = ndgrid(rng{ones(d, 1)});
    grid = cell2mat(permute(grid, [2:d+1 1]));
    grid = reshape(grid, [N^d d]);

    f = exp(2*pi*i*grid*omega/N)*alpha(:);

    f = reshape(f, [N*ones(1, d) 1]);
end
