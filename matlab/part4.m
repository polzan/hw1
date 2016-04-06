close all; clear all; clc;

K = 1000;
N = 2; % Must be the same one of part3
sigma2w = 0.0002;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);

r_estim = autocorrelation_unbiased(x, N-1);
r_theo = theoretical_autocorr(sigma2w, N-1);
R = transpose(toeplitz(r_estim));
lambdas = eig(R);

mu_tilde = 0.175;
mu = mu_tilde/(N*r_estim(1)); %the numerator is mu-tilde
fprintf('Setting mu tilde to %f\n', mu_tilde);
fprintf('mu is %f\n', mu);
fprintf('To be stable it should be under %f (eigenvalues)\n', 2/max(lambdas));
fprintf('To be stable it should be under %f (power)\n', 2/(N*r_estim(1)));

[c, e, y, k_c, k_ey] = lms_predictor(x, N, mu);

fprintf('The value of c at the last iteration is\n');
disp(c(:,K-N+1).');

% Calculate the convergence values for J and c
[a, Jmin] = ar_model(x, N, 'theoretical', sigma2w);
c_opt = -a(2:N+1);
Jinf = 2/(2 - mu*N*r_theo(1)) * Jmin;

figure;
subplot(1,2,1);
hold on;
plot(k_c, real(c(1,:)));
plot(0:K-1, real(c_opt(1)) * ones(K, 1), '--');
ylabel('Real part');
xlabel('k');
xlim([0 150]);
subplot(1,2,2);
hold on;
plot(k_c, imag(c(1,:)));
plot(0:K-1, imag(c_opt(1)) * ones(K, 1), '--');
ylabel('Imaginary part');
xlabel('k');
xlim([0 150]);
print('coeff_lms_c1', '-depsc');

figure;
subplot(1,2,1);
hold on;
plot(k_c, real(c(2,:)));
plot(0:K-1, real(c_opt(2)) * ones(K, 1), '--');
ylabel('Real part');
xlabel('k');
xlim([0 150]);
subplot(1,2,2);
hold on;
plot(k_c, imag(c(2,:)));
plot(0:K-1, imag(c_opt(2)) * ones(K, 1), '--');
ylabel('Imaginary part');
xlabel('k');
xlim([0 150]);
print('coeff_lms_c2', '-depsc');

% figure;
% hold on;
% plot(k, real(x));
% plot(k_ey, real(y));

% simulate to find the MSE
tries=200;
es = zeros(K-N, tries);
for i=1:200
    [x,k] = generate_x(sigma2w, K, 'shuffle');
    [~, e_i] = lms_predictor(x, N, mu);
    es(:,i) = e_i;
end

mse = mean(abs(es).^2, 2);
Jex = Jinf -  Jmin;
fprintf('The excess MSE is %f\n', Jex);
fprintf('The MSD is %f\n', Jex/Jmin);

figure;
hold on;
plot(k_ey, 10*log10(abs(e).^2));
plot(k_ey, 10*log10(mse));
plot(0:K-1, 10*log10(Jmin) * ones(K, 1), '--');
plot(0:K-1, 10*log10(Jinf) * ones(K, 1), '--');
legend('|e|^2', 'E[|e|^2]', 'Jmin', 'Jinf');
ylabel('dB');
xlabel('k');
xlim([0 150]);
print('lms_mse', '-depsc');