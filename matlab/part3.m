close all; clear all; clc;

% params 
K = 1000;
sigma2w = 0.0002;

load 'rng_seed'
[x, k] = generate_x(sigma2w, K, rng_seed);
Fs = 1;
f = linspace(0, Fs, K);

%predictor
N=3;
[c, a] = predictor(x,N,sigma2w);

fprintf('The filter coefficents are:\n');
disp(c);

figure;
[z, p, k] = tf2zpk(a, 1);
zplane(z,p);
print('predictor_zp','-depsc')
