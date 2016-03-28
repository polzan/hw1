close all; clear all; clc;

K = 10000;
N = 4; % Must be the same one of part3
sigma2w = 0.0002;
x = generate_x(sigma2w, K); % Must be the same realization of part3

mu = 0.000005; % ? Max 2/M = 10000 ?
% Skip the first N-1 samples: we need a full vector x_k
c = zeros(N, K-N+1);
y = zeros(K-N+1, 1);
e = zeros(K-N+1, 1);
for k=1:K-N % -1 because we need the last d(k)
    x_k = x(k:k+N-1);
    y(k) = x_k.'*c(:,k);
    e(k) = x(k+N) - y(k); % should predict the next sample
    c(:,k+1) = c(:,k) + mu*e(k)*conj(x_k);
end

fprintf('The value of c at the last iteration is\n');
disp(c(:,K-N+1).');

figure;
plot(1:K-N+1, real(c(1,:)));
