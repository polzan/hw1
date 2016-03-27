function [r, n] = theoretical_autocorr(sigma2w, N)
% fixed params
f_1 = 0.125;
f_2 = 0.8;

n = transpose(0:N-1);
r = exp(1j*2*pi*f_1*n) + exp(1j*2*pi*f_2*n);
r(1) = r(1) + 2*sigma2w;
end
