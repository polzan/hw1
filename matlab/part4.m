close all; clear all; clc;

K = 1000;
N = 2; % Must be the same one of part3
sigma2w = 0.0002;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);

r_theo = theoretical_autocorr(sigma2w, N-1);
R = transpose(toeplitz(r_theo));
lambdas = eig(R);

%mu = 2/(max(lambdas) + min(lambdas)); % ? Max 2/M = 10000 ?
mu = 0.5/(N*r_theo(1));     %the numerator is mu-tilde the quantity we have to report in the solution
fprintf('Setting mu to %f\n', mu);
%fprintf('To be stable it should be under %f\n', 2/max(lambdas));
fprintf('To be stable it should be under %f\n', 2/(N*r_theo(1)));

[c, e, y, k_c, k_ey] = lms_predictor(x, N, mu);

fprintf('The value of c at the last iteration is\n');
disp(c(:,K-N+1).');

figure;
subplot(1,2,1);
plot(k_c, real(c(1,:)));
xlim([0 200]);
subplot(1,2,2);
plot(k_c, imag(c(1,:)));
xlim([0 200]);

figure;
subplot(1,2,1);
plot(k_c, real(c(2,:)));
xlim([0 200]);
subplot(1,2,2);
plot(k_c, imag(c(2,:)));
xlim([0 200]);


% figure;
% hold on;
% plot(0:K-1, real(x));
% plot(N:K-1, real(y));

tries=200;
es = zeros(K-N, tries);
for i=1:200
    [x,k] = generate_x(sigma2w, K, 'shuffle');
    [~, e_i] = lms_predictor(x, N, mu);
    es(:,i) = e_i;
end

mse = mean(abs(es).^2, 2);

figure;
hold on;
plot(k_ey, 10*log10(abs(e).^2));
plot(k_ey, 10*log10(mse));
