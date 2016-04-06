close all; clear all; clc;

K = 1000;
sigma2w = 0.0002;

load 'rng_seed_part5'
[x, k] = generate_x(sigma2w, K, rng_seed_part5);

N = 64;
circle_thresh = 1e-2;

r = autocorrelation_unbiased(x, 0);
mu_tilde = 1;
mu = mu_tilde/(N*r(1)); %the numerator is mu-tilde
fprintf('Setting mu tilde to %f\n', mu_tilde);
fprintf('mu is %f\n', mu);

[c, e, y, k_c, k_ey] = lms_predictor(x, N, mu);

a = [1; -c(:,length(k_c))];
[z,p,~] = tf2zpk(1,a);

figure;
zplane(z,p);
print('part5_poles', '-depsc');

i_poles = find(abs(abs(p) - 1) <= circle_thresh);
poles = p(i_poles);

thetas = angle(poles);
i_less = find(thetas < 0);
thetas(i_less) = thetas(i_less) + 2*pi;
freqs = thetas / (2*pi);

fprintf('The estimated frequencies are:\n');
disp(freqs);
