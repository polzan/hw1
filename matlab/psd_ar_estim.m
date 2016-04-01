function psd = psd_ar_estim(x, N, f, Fs)
[a_tf, sigma2w_estim] = ar_model(x, N);
H = freqz(1, a_tf, 'whole', f, Fs);
psd = sigma2w_estim .* abs(H).^2;
end
