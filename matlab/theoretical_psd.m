function [psd, f] = theoretical_psd(sigma2w, f, Fs)
% params
f_1 = 0.125;
f_2 = 0.8;

K = length(f);
F = 1/K;
psd = 2 * sigma2w * ones(K, 1); % Power of the complex gaussian
i_1 = round(f_1 / F);
i_2 = round(f_2 / F);
psd(i_1) = psd(i_1) + K; % delta as a triangle 2F by K
psd(i_2) = psd(i_2) + K;
end
