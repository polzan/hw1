function r = autocorrelation_biased(x, N, whole)
if nargin < 2
    r_u = autocorrelation_unbiased(x);
    N = length(r_u) -1;
elseif nargin < 3
    r_u = autocorrelation_unbiased(x, N);
    %whole = false;
else
    r_u = autocorrelation_unbiased(x, N, whole);
end
K = length(x);
n = transpose(0:N);
r = (1 - n/K) .* r_u;
end
