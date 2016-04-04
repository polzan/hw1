close all; clear all; clc;

% params 
K = 1000;
sigma2w = 1.26;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);
Fs = 1;
f = linspace(0, Fs, K);

psd_theor = theoretical_psd(sigma2w, f, Fs);

% Periodogram
psd_per = 1/K * abs(fft(x)).^2;

% Welch
D = 100;
S = 25;
win_welch = rectwin(D);
psd_welch = psd_welch_estim(x, win_welch, D, S, f, Fs);

% Correlogram
r = autocorrelation_unbiased(x);
win_corr = blackman(K);
psd_corr = 1/Fs * fft(win_corr .* r);

% AR model
N = 4;
psd_ar = psd_ar_estim(x, N, f, Fs);

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
ylim([-10 40]);
xlim([0 1]);
