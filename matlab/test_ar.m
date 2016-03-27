clear all;
close all;

N = 16;
sigma2w = 1;

r_whole = theoretical_autocorr(sigma2w, N+1);
r_trunc = r_whole(1:N);
R = transpose(toeplitz(r_trunc));

assert(all(R(:,1) == r_trunc));
assert(all(R(:,2) == [conj(r_trunc(2)); r_trunc(1:N-1)]));
assert(all(R(:,N) == flip(conj(r_trunc))));

% fprintf('The determinant of R is %f\n', det(R));
r = r_whole(2:N+1);
%a = -inv(R) * r;
a = linsolve(R, -r);



sigma2w_estim = r_whole(1) + (r'*a);
fprintf('The AR model estimates a noise power of %f\n', sigma2w_estim);
sigma2w_estim = abs(real(sigma2w_estim));

a_whole = [1; a];
f = linspace(0, 1, 1000);
H = freqz(1, a_whole, 'whole', f, 1);

psd = sigma2w_estim .* abs(H).^2;

% = linspace(0, 1, 1000);
plot(f, 10*log10(psd));