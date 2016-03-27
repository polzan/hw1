function psd = psd_ar_estim(x, N, f, Fs)
r_whole = autocorrelation_unbiased(x);
r_trunc = r_whole(1:N);
R = transpose(toeplitz(r_trunc));

assert(all(R(:,1) == r_trunc));
assert(all(R(:,2) == [conj(r_trunc(2)); r_trunc(1:N-1)]));
assert(all(R(:,N) == flip(conj(r_trunc))));

r = r_whole(2:N+1);
%a = -inv(R) * r;
linsolve_opts.SYM = true;
a = linsolve(R, -r, linsolve_opts);

sigma2w_estim = r_whole(1) + (r'*a);
fprintf('The AR model estimates a noise power of %f\n', sigma2w_estim);
sigma2w_estim = abs(real(sigma2w_estim));

a_tf = [1; a];
H = freqz(1, a_tf, 'whole', f, Fs);
psd = sigma2w_estim .* abs(H).^2;
end
