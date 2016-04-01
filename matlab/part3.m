close all; clear all; clc;

% params 
K = 1000;
sigma2w = 0.0002;


load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);
Fs = 1;
f = linspace(0, Fs, K);
% Theoretical psd: F[r(n) = e^j... + e^j... + sigma2*delta] = delta(f - f_1)
% + delta(f - f_2) + sigma2w
psd_theor = 2 * sigma2w * ones(K, 1);
psd_theor(round(f_1 * K)) = 2 * sigma2w + K;
psd_theor(round(f_2 * K)) = 2 * sigma2w + K;

% Periodogram
X = fft(x); % Tc = 1
psd_per = 1/K * abs(X).^2;
%psd_per_matlab = periodogram(x, rectwin(K), K).';

% Welch
D = 70;
S = 20;
win_welch = blackman(D);
psd_welch = psd_welch_estim(x, win_welch, D, S);

% Correlogram
r = autocorrelation_unbiased(x);
win_corr = blackman(K);
psd_corr = fft(win_corr .* r);

% AR model
N = 16;
psd_ar = psd_ar_estim(x, N);
%psd_ar_matlab = pcov(x, N, K);

%predictor
N=4;
[c, fn, a] = predictor(x,N);
z = tf2zpk(a,1);
zplane(z)   %something wrong (not a minimum phase filter)


figure;
hold all;
plot(f, 10*log10(abs(psd_corr)), 'LineWidth', 1.5);
plot(f, 10*log10(abs(psd_per)), 'LineWidth', 1.5);
plot(f, 10*log10(abs(psd_welch)), 'LineWidth', 1.5);
plot(f, 10*log10(abs(psd_ar)), 'LineWidth', 1.5);
plot(f, 10*log10(abs(psd_theor)), 'LineWidth', 1.5);

legend('Correlogram', 'Periodogram', 'Welch', 'AR', 'Theoretical');
ylabel('dB');
xlabel('Hz');
%ylim([0 50]);