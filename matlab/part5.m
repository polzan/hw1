close all; clear all; clc;

K = 1000;
N = 20; 
sigma2w = 0.0002;

load 'rng_seed_part5'
[x, k] = generate_x(sigma2w, K, rng_seed_part5);

D = 5; % Everything > 0 should be fine: the broadband sig. is white
x_delayed = x(D+1:K);
k_delayed = D:K-1;

r = autocorrelation_unbiased(x_delayed, 0);

mu_tilde = 0.25;
mu = mu_tilde/(N*r(1)); %the numerator is mu-tilde
fprintf('Setting mu tilde to %f\n', mu_tilde);
fprintf('Setting mu to %f\n', mu);

[c, e, y, k_c, k_ey] = lms_filter(x(1:K-D+1), x_delayed, N, mu);

% figure;
% plot(k_ey, abs(autocorrelation_unbiased(y, length(k_ey)-1)));
% 
% figure;
% plot(k_ey, abs(autocorrelation_unbiased(e, length(k_ey)-1)));

Nar = N;
[a, ~] = ar_model(y, Nar);
[z,p,~] = tf2zpk(1,a);

figure;
zplane(z,p);
print('part5poles', '-depsc');

i_poles = find(abs(abs(p) - 1) <= 1e-3);
poles = p(i_poles);

thetas = angle(poles);
i_less = find(thetas < 0);
thetas(i_less) = thetas(i_less) + 2*pi;
freqs = thetas / (2*pi);

fprintf('The estimated frequencies are:\n');
disp(freqs);
