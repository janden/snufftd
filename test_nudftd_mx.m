% Script to compare the MATLAB and MEX implementations of the NUDFT.

% Parameters.
N = 64;
n = 2^10;
d = 3;

% Generate sample data.
rand('state', 0);
omega = N/2*(rand(n, d)-0.5);
alpha = rand(n, 1) + 1i*rand(n, 1);

% MATLAB implementation: nudftd.
tmr = tic;
f0 = nudftd(N, omega', alpha);
tm0 = toc(tmr);

% MEX implementation in C: nudftd_mx.
tmr = tic;
f1 = nudftd_mx(N, omega, alpha);
tm1 = toc(tmr);

% Output results.
fprintf('%-20s %f s\n', 'nudftd:', tm0);
fprintf('%-20s %f s\n', 'nudftd_mx:', tm1);
fprintf('%-20s %g\n', 'Relative error:', norm(f0(:)-f1(:))/norm(f0(:)));
