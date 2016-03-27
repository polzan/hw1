close all; clear all; clc;

% params 
K = 1000;
sigma2w = 0.0002;
f_1 = 0.125;
f_2 = 0.8;

% r.v.s
phi_1 = 2*pi*rand(1);
phi_2 = 2*pi*rand(1);
w = sqrt(sigma2w) * (randn(K, 1) + 1j * randn(K, 1));

% generate the r.p. x
k = transpose(0:(K - 1));
x = exp(1j*2*pi*f_1*k + phi_1) + exp(1j*2*pi*f_2*k + phi_2) + w;
%x = exp(1j*2*pi*f_1*k);

N = 30;

r = autocorrelation_unbiased(x);
R_N = transpose(toeplitz(r(1:N)));
r_N = r(2:N+1);

linsolve_opts.SYM = true;
c_opt = linsolve(R_N, r_N, linsolve_opts);

a_N = [1; -c_opt];

zplane(1,a_N);

f_N = filter(1, a_N, x);

figure;
subplot(2,1,1);
plot(0:K-1, 20*log10(abs(fft(f_N))));
subplot(2,1,2);
plot(0:K-1,20*log10(abs(fft(x))));

