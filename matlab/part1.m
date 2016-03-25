close all;
clear all;

% params 
K = 1000;
sigma2w = 1.26;
f_1 = 0.125;
f_2 = 0.8;

% r.v.s
phi_1 = 2*pi*rand(1);
phi_2 = 2*pi*rand(1);
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
D = 70;
S = 20;
win_welch = blackman(D);
Mw = 1/D * (win_welch'*win_welch);
%psd_welch_matlab = pwelch(x, win_welch, S, K);
N_s = floor((K - D) / (D - S) + 1);
psd_w = zeros(K, N_s);
for s = 0:N_s - 1
    x_s = win_welch .* x((s*(D-S) + 1):(s*(D-S) + D));
    X_s = fft(x_s, K);
    psd_w(:,s+1) = 1/(D*Mw) * abs(X_s).^2;
end

psd_welch = mean(psd_w, 2);

% Correlogram
r = autocorrelation_unbiased(x);
win_corr = blackman(K);
psd_corr = fft(win_corr .* r);

% AR model
N = 25;
r_whole = autocorrelation_unbiased(x);
r_trunc = r(1:N);
R = toeplitz(r_trunc).';

assert(all(R(:,1) == r_trunc));
assert(all(R(:,2) == [conj(r_trunc(2)); r_trunc(1:N-1)]));
assert(all(R(:,N) == flip(conj(r_trunc))));

fprintf('The determinant of R is %f\n', det(R));
r = r_whole(2:N+1);
a = -inv(R) * r;

sigma2w_estim = r_whole(1) + (r'*a);
fprintf('The AR model estimates a noise power of %f\n', sigma2w_estim);

H = freqz(1, a, 1000, 1, 'whole');

%figure;
%plot(f, 20*log10(abs(H)));

psd_ar = 2*sigma2w .* abs(H).^2;

% x2 = filter(1, a, w);
% figure;
% plot(k, 20*log10(abs(x2 - x)));

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
ylim([0 50]);
