function psd = psd_ar_estim(x, N)
r_whole = autocorrelation_unbiased(x);
r_trunc = r_whole(1:N);
R = toeplitz(r_trunc).';

% assert(all(R(:,1) == r_trunc));
% assert(all(R(:,2) == [conj(r_trunc(2)); r_trunc(1:N-1)]));
% assert(all(R(:,N) == flip(conj(r_trunc))));

% fprintf('The determinant of R is %f\n', det(R));
r = r_whole(2:N+1);
a = -inv(R) * r;

sigma2w_estim = r_whole(1) + (r'*a);
%fprintf('The AR model estimates a noise power of %f\n', sigma2w_estim);
sigma2w_estim = abs(real(sigma2w_estim));

H = freqz(1, a, 1000, 1, 'whole');

psd = sigma2w_estim .* abs(H).^2;
end
