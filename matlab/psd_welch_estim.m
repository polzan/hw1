function psd = psd_welch_estim(x, w, D, S, f, Fs)
K = length(x);
Ndft = length(f);
Mw = 1/D * (w'*w);
N_s = floor((K - D) / (D - S) + 1);
psd_chunk = zeros(Ndft, N_s);
for s = 0:N_s - 1
    x_s = w .* x((s*(D-S) + 1):(s*(D-S) + D));
    psd_chunk(:,s+1) = 1/(D*Mw) * abs(fft(x_s, Ndft)).^2;
end
psd = mean(psd_chunk, 2);
end
