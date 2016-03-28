close all; clear all; clc;

[x, k, w] = generate_x(0.0002, 1000);

N = 16;
K = length(x);
%r_whole=theoretical_autocorr(0.0002, N+1);
r_whole = autocorrelation_unbiased(x);
r_n = r_whole(2:N+1);
R_n = transpose(toeplitz(r_whole(1:N))); 
c_opt = inv(R_n)*r_n;
a = [1;-c_opt];
f_n = conv(x,a, 'same');

figure;
[H, f] = freqz(1, a, 'whole', 1000, 1);
subplot(2,1,1);
plot(f, 20*log10(abs(H)));
subplot(2,1,2);
plot(f, angle(H) / pi);
figure;
[z, p] = tf2zpk(1, a);
zplane(1,a);

figure;
diff = abs(w - f_n);
plot(k, 20*log10(diff));
