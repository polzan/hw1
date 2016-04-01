function r = autocorrelation_biased(x, N)
if nargin < 2
    r_u = autocorrelation_unbiased(x);
    N = length(r_u) -1;
else
r_u = autocorrelation_unbiased(x, N);
end
K = length(x);
n = transpose(0:N);
r = (1 - n/K) .* r_u;
end
