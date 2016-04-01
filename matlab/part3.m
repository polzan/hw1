close all; clear all; clc;

% params 
K = 1000;
sigma2w = 0.0002;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);
Fs = 1;
f = linspace(0, Fs, K);

%predictor
N=2;
[c, fn, a] = predictor(x,N);

figure;
[H, f] = freqz(1, a, 'whole', 1000, 1);
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
subplot(2,1,2);
plot(f, angle(H) / pi);

figure;
%[z, p, k] = tf2zpk(1, a);
zplane(1,a);
