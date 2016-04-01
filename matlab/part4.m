close all; clear all; clc;

K = 1000;
N = 4; % Must be the same one of part3
sigma2w = 0.0002;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);

r = theoretical_autocorr(sigma2w, N);
R = transpose(toeplitz(r));
lambda = eig(R);

%mu = 2/(max(lambda) + min(lambda)); % ? Max 2/M = 10000 ?
%mu = 1e-5;
mu = 0.25;
% Skip thefirst N-1 samples: we need a full vector x_k
c = zeros(N, K-N+1);
y = zeros(K-N+1, 1);
e = zeros(K-N+1, 1);
for k=1:K-N % -1 because we need the last d(k)
    x_k = x(k : k + N - 1);
    y(k) = x_k.'*c(:,k);
    e(k) = x(k+N) - y(k); % should predict the next sample
    c(:,k+1) = c(:,k) + mu*e(k) .* conj(x_k);
    %fprintf('deltac:\n');
    %disp(mu*e(k) .* conj(x_k));
    %pause;
    %c(:,k+1) = [1; zeros(N-1,1)];
end

fprintf('The value of c at the last iteration is\n');
disp(c(:,K-N+1).');

figure;
plot(1:K-N+1, real(c(1,:)));

figure;
hold on;
plot(0:K-1, real(x));
plot(N-1:K-1, real(y));

r = theoretical_autocorr(sigma2w, N);
R = transpose(toeplitz(r));
