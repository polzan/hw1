close all; clear all; clc;

K = 1000;
N = 2; % Must be the same one of part3
sigma2w = 0.0002;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);

D = 5;
x_delayed = x(D+1:K);
k_delayed = D:K-1;

r = theoretical_autocorr(sigma2w, N-1);
R = transpose(toeplitz(r));
lambdas = eig(R);

mu_tilde = 0.25;
mu = mu_tilde/(N*r(1)); %the numerator is mu-tilde
fprintf('Setting mu tilde to %f\n', mu_tilde);
fprintf('Setting mu to %f\n', mu);

[c, e, y, k_c, k_ey] = lms_filter(x(1:K-D+1), x_delayed, N, mu);

figure;
plot(0:K-D-2, abs(autocorrelation_unbiased(y, K-D-2)));

figure;
plot(0:K-D-2, abs(autocorrelation_unbiased(e, K-D-2)));

Nar = 6;
[a, sigma2w] = ar_model(y, Nar);
[z,p,~] = tf2zpk(1,a);

figure;
zplane(z,p);