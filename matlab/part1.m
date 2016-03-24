close all;
clear all;

% params 
K = 1000;
sigma2w = 1.26;
f_1 = 0.125;
f_2 = 0.8;

% r.v.s
phi_1 = 2*pi*rand(K, 1);
phi_2 = 2*pi*rand(K, 1);
w = sqrt(sigma2w) * (randn(K, 1) + 1j * randn(K, 1));

% generate the r.p. x
k = transpose(0:(K - 1));
x = exp(1j*2*pi*f_1*k + phi_1) + exp(1j*2*pi*f_2*k + phi_2) + w;

f = linspace(0, 1, K);
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
D = 50;
S = 10;
win = rectwin(D);
Mw = 1/D * (win'*win);
%psd_welch_matlab = pwelch(x, w, S, K);
N_s = floor((K - D) / (D - S) + 1);
psd_w = zeros(K, N_s);
for s = 0:N_s - 1
    x_s = win .* x((s*(D-S) + 1):(s*(D-S) + D));
    X_s = fft(x_s, K);
    psd_w(:,s+1) = 1/(D*Mw) * abs(X_s).^2;
end

psd_welch = mean(psd_w, 2);

% Correlogram
r = autocorrelation_unbiased(x);
win_corr = rectwin(K);
psd_corr = fft(win_corr .* r);

% AR model


% figure;
% hold all;
% plot(f, 10*log10(psd_w(:,1)));
% plot(f, 10*log10(psd_w(:,2)));
% plot(f, 10*log10(psd_w(:,3)));

figure;
hold all;
plot(f, 10*log10(psd_theor));
plot(f, 10*log10(psd_per));
%plot(f, 10*log10(psd_per_matlab));
%plot(f, 10*log10(psd_welch_matlab));
plot(f, 10*log10(psd_welch));
plot(f, 10*log10(psd_corr));


