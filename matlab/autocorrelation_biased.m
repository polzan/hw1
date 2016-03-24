function r = autocorrelation_biased(x)
r_u = autocorrelation_unbiased(x);
K = length(x);
n = transpose(0:K-1);
r = (1 - abs(n)/K) .* r_u;
end
