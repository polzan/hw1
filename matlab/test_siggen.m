close all; clear all; clc;

seed='shuffle';
K = 1000;
sigma2w = 1.26;

% params
f_1 = 0.125;
f_2 = 0.8;

% r.v.s
rngstate = rng(seed);
phi_1 = 2*pi*rand(1);
phi_2 = 2*pi*rand(1);
w = sqrt(sigma2w/2) * (randn(K, 1) + 1j * randn(K, 1));  %sigma2w_i = sigma2w/2 (complex gaussian)
rng(rngstate);

assert(phi_1 <= 2*pi);
assert(phi_1 >= 0);
assert(phi_2 <= 2*pi);
assert(phi_2 >= 0);

figure;
subplot(1,2,1);
histogram(real(w));
subplot(1,2,2);
histogram(imag(w));

% generate the r.p. x
k = transpose(0:(K - 1));
exp1 = exp(1j*(2*pi*f_1*k + phi_1));
exp2 = exp(1j*(2*pi*f_2*k + phi_2));

figure;
subplot(1,2,1);
plot(real(exp1), imag(exp1));
subplot(1,2,2);
plot(real(exp2), imag(exp2));

x = exp1 + exp2 + w;

figure;
scatter(real(x), imag(x));
