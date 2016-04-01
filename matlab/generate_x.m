function [x, k, w, phi_1, phi_2] = generate_x(sigma2w, K, seed)
if nargin < 3
    seed = 'shuffle';
end

% params
f_1 = 0.125;
f_2 = 0.8;

% r.v.s
rngstate = rng(seed);
phi_1 = 2*pi*rand(1);
phi_2 = 2*pi*rand(1);
w = sqrt(sigma2w/2) * (randn(K, 1) + 1j * randn(K, 1));  %sigma2w_i = sigma2w/2 (complex gaussian)
rng(rngstate);

% generate the r.p. x
k = transpose(0:(K - 1));
x = exp(1j*(2*pi*f_1*k + phi_1)) + exp(1j*(2*pi*f_2*k + phi_2)) + w;
end
