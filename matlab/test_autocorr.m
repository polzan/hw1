close all; clear all; clc;

sigma2w = 1.26;

[x, k] = generate_x(sigma2w, 10001);
r = theoretical_autocorr(sigma2w, 10001);
r_estim = autocorrelation_biased(x);
%r_estim = r_estim(1:10000);

r_matlab = autocorr(x, 10000);

figure;
hold all;
plot(k, abs(r)./abs(r(1)));
plot(k, abs(r_estim)./abs(r_estim(1)));
plot(k, abs(r_matlab)./abs(r_matlab(1)));

iter = 100;
K = 10000;
r = zeros(K, iter);
for i=1:iter
   x = generate_x(1.5, K);
   r(:,i) = autocorrelation_biased(x);
   %r(:,i) = autocorr(x, K-1);
end

r_mean = mean(r, 2);

figure;
plot(0:K-1, abs(r_mean));
