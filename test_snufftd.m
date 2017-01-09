% Test script for d-dimensional NUFFT and S-NUFFT.

% Set parameters.
N = 8;
n = 5;
d = 3;

% Generate data.
rand('state', 0);

omega = N*(rand(n, d)-0.5);
alpha = rand(n, 1) + 1i*rand(n, 1);

% Non-uniform DFT.
tmr = tic;
f0 = nudftd(N, omega, alpha);
tm0 = toc(tmr);
fprintf('%-15sTime: %15f s\n', 'NUDFT', tm0);

% Standard non-uniform FFT.
tmr = tic;
f1 = nufftd(N, omega, alpha);
tm1 = toc(tmr);
err1 = norm(f0(:)-f1(:));
fprintf('%-15sTime: %15f s    Error: %15g\n', 'NUFFT', tm1, err1);

% MEX-augmented non-uniform FFT.
tmr = tic;
f2 = nufftd(N, omega, alpha, [], [], [], true);
tm2 = toc(tmr);
err2 = norm(f0(:)-f2(:));
fprintf('%-15sTime: %15f s    Error: %15g\n', 'NUFFT (MEX)', tm2, err2);

% Compact non-uniform FFT.
tmr = tic;
f3 = snufftd(N, omega, alpha);
tm3 = toc(tmr);
err3 = norm(f0(:)-f3(:));
fprintf('%-15sTime: %15f s    Error: %15g\n', 'CNUFFT', tm3, err3);

if exist('nufft1d1') && d <= 3
    if d == 1
        tmr = tic;
        f4 = n*nufft1d1(n, ...
            2*pi/N*omega(:,1), ...
            alpha, 1, 1e-10, N);
        tm4 = toc(tmr);
        fun4 = 'nufft1d1';
    elseif d == 2
        tmr = tic;
        f4 = n*nufft2d1(n, ...
            2*pi/N*omega(:,1), 2*pi/N*omega(:,2), ...
            alpha, 1, 1e-10, N, N);
        tm4 = toc(tmr);
        fun4 = 'nufft2d1';
    elseif d == 3
        tmr = tic;
        f4 = n*nufft3d1(n, ...
            2*pi/N*omega(:,1), 2*pi/N*omega(:,2), 2*pi/N*omega(:,3), ...
            alpha, 1, 1e-10, N, N, N);
        tm4 = toc(tmr);
        fun4 = 'nufft3d1';
    end
    err4 = norm(f0(:)-f4(:));
    fprintf('%-15sTime: %15f s    Error: %15g\n', fun4, tm4, err4);
end
