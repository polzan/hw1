function psd = psd_welch_estim(x, w, D, S)
K = length(x);
Mw = 1/D * (w'*w);
N_s = floor((K - D) / (D - S) + 1);
psd_w = zeros(K, N_s);
for s = 0:N_s - 1
    x_s = w .* x((s*(D-S) + 1):(s*(D-S) + D));
    X_s = fft(x_s, K);
    psd_w(:,s+1) = 1/(D*Mw) * abs(X_s).^2;
end
psd = mean(psd_w, 2);
end

